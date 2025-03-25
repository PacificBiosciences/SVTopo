use std::{
    collections::{BTreeMap, HashMap, HashSet},
    time::SystemTime,
};

use log::debug;

use crate::{
    containers::{
        ComplexSVBlock, ComplexSVCalls, Coordinate, EventGraph, FwdStrandSplitReadSegment,
        Orientation,
    },
    utils::MAX_GRAPH_SIZE,
};

/// Starting with event graphs of connected breakends, generate annotated version with
/// unambiguously connected breakends. Returns a vector of ComplexSvBlock vectors,
/// where there is one vector entry for each complex sv graph. Each entry contains a vector of coverage
/// data in ComplexSvBlocks, each representing the data for either one breakend or two connected breakends
pub fn annotate_graphs(
    event_graphs: &Vec<EventGraph>,
    clip_coordinates: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
) -> ComplexSVCalls {
    let start_time = SystemTime::now();
    let mut annotated_event_graphs = Vec::new();

    for event_graph in event_graphs {
        if event_graph.graph.len() > MAX_GRAPH_SIZE {
            log::debug!(
                "Event graph too large, skipped. Total number of blocks in this event: {}.",
                event_graph.graph.len(),
            );
            let ends = get_event_ends(&event_graph.graph);
            if let (Some(start), Some(end)) = ends {
                log::debug!("Skipped event: {} -> {}", start, end,);
            }
            let mut keys: Vec<&u32> = event_graph.graph.keys().collect();
            keys.sort();
            let start_idx = keys.first().unwrap();
            let end_idx = keys.last().unwrap();
            let starts = event_graph.graph.get(start_idx).unwrap();
            let ends = event_graph.graph.get(end_idx).unwrap();
            let start = starts.first().unwrap();
            let end = ends.first().unwrap();

            log::debug!(
                "Event graph too large, skipped. Total number of blocks in this event: {}. From {} to {}",
                event_graph.graph.len(),
                start,
                end,
            );
            continue;
        }
        let mut annotated_event_graph: Vec<ComplexSVBlock> =
            annotate_graph(&event_graph.graph, clip_coordinates);

        // Last block won't be spanned if an unspanned block is last, which can
        // happen in cases where a phasing failure splits up blocks incorrectly.
        // It will also be unspanned if there is no last block, aka an empty event graph.
        let mut last_block_is_spanned = false;
        if let Some(last_block) = annotated_event_graph.last() {
            last_block_is_spanned = !last_block.coverages.is_empty();
        }
        if last_block_is_spanned {
            if let Some(first_block) = annotated_event_graph.first() {
                if first_block.coverages.is_empty() {
                    continue;
                }
            }
            add_orientations(&mut annotated_event_graph);
            annotated_event_graphs.push(annotated_event_graph);
        }
    }
    debug!("Annotation: {}s", start_time.elapsed().unwrap().as_secs());
    debug!("{} annotated blocks", annotated_event_graphs.len());
    ComplexSVCalls::new(annotated_event_graphs)
}

/// Utility primarily for logging, gets the first and last blocks from an event graph
fn get_event_ends(
    event_graph: &HashMap<u32, Vec<Coordinate>>,
) -> (Option<&Coordinate>, Option<&Coordinate>) {
    let mut coords = Vec::new();
    for coordinate_level in event_graph.values() {
        for coord in coordinate_level.iter() {
            coords.push(coord);
        }
    }
    coords.sort();
    (coords.first().copied(), coords.last().copied())
}

/// Annotate a single complex event graph with coverages and updated order indices,
/// combining those that can be joined by alignments.
/// Update the sample order indices when combined.
///
/// **Basic steps are**:
/// 1) Loop through the coordinates at each sample order index in the event graph.
/// 2) Add start: Find coverage going backward (upstream) from the first block and add it as the start.
/// 3) Add unspanned: Find coverage for unspanned connections going back from the current coordinate to
///     any at the previous level. Add these as unspanned connections.
/// 4) Add spanned: Find coverage for spanned connections going forward from the current coordinate.
/// 5) Add tied connections: Find coverage for unspanned connections between this coordinate
///     and any others at the same sample order index, then add them as unspanned connections.
/// 6) When all levels have been traversed, add a final spanned coordinate from the last one going downstream.
fn annotate_graph(
    event_graph: &HashMap<u32, Vec<Coordinate>>,
    clip_coordinates: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
) -> Vec<ComplexSVBlock> {
    let mut annotated_event_graph: Vec<ComplexSVBlock> = Vec::new();
    if let Some(max_idx) = event_graph.keys().max() {
        // the complex event start is in forward orientation. Each directly connected event orientation after
        // can be inferred from it, but when the chain breaks the orientation becomes unknowable
        let max_idx = *max_idx as i64;
        let mut already_processed_coordinates: HashSet<Coordinate> = HashSet::new();

        let mut last_event_spanned = false;
        let mut order_scaler = 0;
        for sample_idx in 0..max_idx + 1 {
            let is_first = sample_idx == 0;
            if let Some(coordinates) = event_graph.get(&(sample_idx as u32)).cloned() {
                let mut prev_coords = Vec::new();
                if !is_first {
                    if let Some(found_prev_coords) = event_graph.get(&(sample_idx as u32 - 1)) {
                        prev_coords.append(&mut found_prev_coords.clone());
                    }
                }
                (last_event_spanned, order_scaler) = add_blocks(
                    &coordinates,
                    event_graph,
                    &mut already_processed_coordinates,
                    is_first,
                    sample_idx == max_idx,
                    clip_coordinates,
                    &mut annotated_event_graph,
                    prev_coords,
                    last_event_spanned,
                    order_scaler,
                    sample_idx,
                );
            }
        }
        add_last_block(
            clip_coordinates,
            &mut annotated_event_graph,
            event_graph,
            max_idx as u32,
        );
    }
    annotated_event_graph
}

/// Loop through and add coordinate blocks that are tied at the same
/// level of the sample order to the annotated event graph.
fn add_blocks(
    coordinates: &[Coordinate],
    event_graph: &HashMap<u32, Vec<Coordinate>>,
    already_processed_coordinates: &mut HashSet<Coordinate>,
    is_first: bool,
    is_last: bool,
    clip_coordinates: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    annotated_event_graph: &mut Vec<ComplexSVBlock>,
    mut prev_coords: Vec<Coordinate>,
    mut last_event_spanned: bool,
    order_scaler: u32,
    sample_idx: i64,
) -> (bool, u32) {
    let mut max_order_scaler = 0;
    for (i, coordinate) in coordinates.iter().enumerate() {
        let mut curr_order_scaler = order_scaler;
        if already_processed_coordinates.contains(coordinate) {
            continue;
        }
        if is_first {
            add_first_block(
                clip_coordinates,
                coordinate,
                annotated_event_graph,
                &mut prev_coords,
            );
            last_event_spanned = true;
        } else {
            let mut unspanned_added = false;
            for (mut unspanned_block, _) in
                get_unspanned_connections(coordinate, &prev_coords, clip_coordinates)
            {
                unspanned_block.sample_order_index = curr_order_scaler + sample_idx as u32;
                annotated_event_graph.push(unspanned_block);
                unspanned_added = true;
                already_processed_coordinates.insert(coordinate.clone());
                last_event_spanned = false;
            }
            let blocks_to_add: Vec<(ComplexSVBlock, Coordinate)> = generate_spanned_blocks(
                sample_idx,
                event_graph,
                clip_coordinates,
                coordinate,
                already_processed_coordinates,
                is_last,
            );
            if last_event_spanned {
                curr_order_scaler += 1;
            }
            for (mut spanned_block, _) in blocks_to_add {
                spanned_block.sample_order_index = curr_order_scaler + sample_idx as u32;
                if unspanned_added {
                    spanned_block.sample_order_index += 1;
                }
                annotated_event_graph.push(spanned_block);
                last_event_spanned = true;
            }
        }
        if curr_order_scaler > max_order_scaler {
            max_order_scaler = curr_order_scaler;
        }

        // if there are multiple tied coordinates at this sample_order_index level and this is the first one,
        // check for unspanned connections that join them
        if coordinates.len() > 1 && i == 0 {
            if let Some(tied_coords) = event_graph.get(&(sample_idx as u32)) {
                let tied_unspanned_blocks =
                    get_unspanned_connections(coordinate, tied_coords, clip_coordinates);
                for (mut tied_unspanned_block, _) in tied_unspanned_blocks {
                    tied_unspanned_block.sample_order_index = curr_order_scaler + sample_idx as u32;
                    annotated_event_graph.push(tied_unspanned_block);
                }
            };
        }
    }
    (last_event_spanned, max_order_scaler)
}

/// Add the first block to the graph
fn add_first_block(
    clip_coordinates: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    coordinate: &Coordinate,
    annotated_graph: &mut Vec<ComplexSVBlock>,
    previous_coordinates: &mut Vec<Coordinate>,
) {
    // check for coverage for a block before this one, which won't start with a break
    let coverages = get_coordinate_alignments_as_coverages(
        clip_coordinates,
        coordinate,
        None,
        false,
        true,
        false,
    );
    let spanning_coordinate_opt: Option<Coordinate> = coordinate_from_coverages(&coverages);
    if let Some(mut spanned_connection) = spanning_coordinate_opt.clone() {
        spanned_connection.variant_ids = coordinate.variant_ids.clone();
        annotated_graph.push(ComplexSVBlock::new(spanned_connection, coverages, 0, "+"));
        previous_coordinates.push(coordinate.clone());
    }
}

/// Add the last block to the graph if missing, starting with unspanned connections to previous blocks.
/// This is only necessary if the last block shares a starting breakpoint with another location in the
/// graph, as that will create an ambigous block that won't be created.
///
/// **Steps:**
///
/// 1. For each coordinate from the un-annotated event graph at the final sample index,
///     perform a fuzzy match check to see if it's been added already (match to previous end coord).
///     If not already added, continue with the following steps.
/// 2. Check for any unspanned connections to previous and add them.
/// 3. Check for spanned connections from this putative end coordinate going downstream, add if found.
fn add_last_block(
    clip_coordinates: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    annotated_event_graph: &mut Vec<ComplexSVBlock>,
    event_graph: &HashMap<u32, Vec<Coordinate>>,
    max_sample_idx: u32,
) {
    let mut current_sample_order_index = 0;
    if let Some(prev_block) = annotated_event_graph.last() {
        current_sample_order_index = prev_block.sample_order_index + 1;
    }
    let mut prev_ends = HashSet::new();
    for block in annotated_event_graph.iter() {
        if block.sample_order_index == current_sample_order_index - 1 {
            prev_ends.insert(block.clone());
        }
    }

    // identify and add unspanned connections to previous
    if let Some(last_coords) = event_graph.get(&max_sample_idx) {
        for last_coord in last_coords {
            if already_added(&prev_ends, last_coord) {
                continue;
            }
            let prev_coords =
                get_previous_blocks(current_sample_order_index, annotated_event_graph);

            // Find and add unspanned connections to prev blocks
            if add_last_unspanned_connection(
                current_sample_order_index,
                last_coord,
                prev_coords,
                clip_coordinates,
                annotated_event_graph,
            ) {
                current_sample_order_index += 1;
            }

            // check for coverage for a block after this one, which won't end with a break
            let coverages = get_coordinate_alignments_as_coverages(
                clip_coordinates,
                last_coord,
                None,
                false,
                false,
                true,
            );
            let spanning_coordinate_opt: Option<Coordinate> = coordinate_from_coverages(&coverages);
            if let Some(mut spanned_connection) = spanning_coordinate_opt {
                spanned_connection.variant_ids = last_coord.variant_ids.clone();
                annotated_event_graph.push(ComplexSVBlock::new(
                    spanned_connection.clone(),
                    coverages,
                    current_sample_order_index,
                    "+",
                ));
            }
        }
    }
}

/// Check if the last coordinate was already added to the annotated event graph
fn already_added(prev_ends: &HashSet<ComplexSVBlock>, last_coord: &Coordinate) -> bool {
    let mut already_added = false;
    for prev_final_block in prev_ends.iter() {
        let mut prev_final_block_end = Coordinate::new(
            prev_final_block.region.start_chrom.clone(),
            prev_final_block.region.start,
        );
        if prev_final_block.orientation == Orientation::Forward {
            prev_final_block_end = Coordinate::new(
                prev_final_block.region.end_chrom.clone(),
                prev_final_block.region.end,
            );
        }

        if last_coord.is_within(&prev_final_block_end) {
            already_added = true;
        }
    }
    already_added
}

/// Get the blocks previous to this sample order index in the graph
fn get_previous_blocks(
    current_sample_order_index: u32,
    annotated_event_graph: &mut [ComplexSVBlock],
) -> Vec<Coordinate> {
    let mut prev_coords = Vec::new();
    for prev_block in annotated_event_graph.iter() {
        if prev_block.sample_order_index == (current_sample_order_index - 1) {
            prev_coords.push(Coordinate::new(
                prev_block.region.start_chrom.clone(),
                prev_block.region.start,
            ));
            prev_coords.push(Coordinate::new(
                prev_block.region.end_chrom.clone(),
                prev_block.region.end,
            ));
        }
    }
    prev_coords
}

fn add_last_unspanned_connection(
    current_sample_order_index: u32,
    last_coord: &Coordinate,
    prev_coords: Vec<Coordinate>,
    clip_coordinates: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    annotated_event_graph: &mut Vec<ComplexSVBlock>,
) -> bool {
    let mut unspanned_connection_found = false;
    let unspanned_connections =
        get_unspanned_connections(last_coord, &prev_coords, clip_coordinates);
    for (mut unspanned_connection, mut _connection_coord) in unspanned_connections {
        if (unspanned_connection.region.start - last_coord.start).abs()
            <= last_coord.confidence_interval.0
            && unspanned_connection.orientation == Orientation::Reverse
            || (unspanned_connection.region.end - last_coord.start).abs()
                <= last_coord.confidence_interval.0
                && unspanned_connection.orientation == Orientation::Forward
        {
            unspanned_connection.sample_order_index = current_sample_order_index;
            annotated_event_graph.push(unspanned_connection);
            unspanned_connection_found = true;
        }
    }
    unspanned_connection_found
}

/// Generate a vector of spanned blocks that this given coordinate is part of. This means
/// generally either the block this coordinate already describes or the block made from this
/// one and a second one, joined by phasing.
fn generate_spanned_blocks(
    sample_idx: i64,
    event_graph: &HashMap<u32, Vec<Coordinate>>,
    clip_coordinates: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    coordinate: &Coordinate,
    already_processed_coordinates: &mut HashSet<Coordinate>,
    is_last: bool,
) -> Vec<(ComplexSVBlock, Coordinate)> {
    let mut blocks_to_add: Vec<(ComplexSVBlock, Coordinate)> = Vec::new();
    // if this coordinate is the first part of a fully spanned block, look ahead and add the second
    // part to this block entry, then include the second part in already_processed_coordinates
    let next_coords_opt = event_graph.get(&(sample_idx as u32 + 1));
    let spanned_blocks = get_spanned_connections_to_next(
        coordinate,
        &next_coords_opt,
        clip_coordinates,
        already_processed_coordinates,
    );
    for spanned_block in spanned_blocks.iter() {
        blocks_to_add.push((spanned_block.clone(), coordinate.clone()));
    }

    // if not already added as part of a combined block or first block, add it now solo
    if sample_idx > 0 && spanned_blocks.is_empty() {
        let block_coverages = get_coordinate_alignments_as_coverages(
            clip_coordinates,
            coordinate,
            None,
            false,
            false,
            is_last,
        );
        let block_coordinate_opt = coordinate_from_coverages(&block_coverages);
        if let Some(mut block_coordinate) = block_coordinate_opt {
            block_coordinate.variant_ids = coordinate.variant_ids.clone();
            let block = ComplexSVBlock::new(block_coordinate, block_coverages, 0, "");
            blocks_to_add.push((block, coordinate.clone()));
        }
    }
    if !blocks_to_add.is_empty() {
        already_processed_coordinates.insert(coordinate.clone());
    }
    blocks_to_add
}

/// Given a coordinate and a list of next coordinates in the graph,
/// created ComplexSVBLock connections for any case where the block
/// includes both the original coordinate and one of the next coordinates.
fn get_spanned_connections_to_next(
    current_coord: &Coordinate,
    next_coords_opt: &Option<&Vec<Coordinate>>,
    clip_coordinates: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    already_processed_coordinates: &mut HashSet<Coordinate>,
) -> Vec<ComplexSVBlock> {
    let mut spanned_blocks: Vec<ComplexSVBlock> = Vec::new();
    if let Some(next_coordinates) = next_coords_opt {
        for next_coord in next_coordinates.iter() {
            let next_coord_opt = Some(next_coord.clone());
            let spanned_coverages = get_coordinate_alignments_as_coverages(
                clip_coordinates,
                current_coord,
                next_coord_opt,
                false,
                false,
                false,
            );
            if !spanned_coverages.is_empty() {
                let spanning_coordinate_opt: Option<Coordinate> =
                    coordinate_from_coverages(&spanned_coverages);
                if let Some(mut spanned_coordinate) = spanning_coordinate_opt {
                    let variant_ids: Vec<String> = current_coord
                        .variant_ids
                        .iter()
                        .chain(next_coord.variant_ids.iter())
                        .cloned()
                        .collect();
                    spanned_coordinate.variant_ids = variant_ids;

                    let new_block =
                        ComplexSVBlock::new(spanned_coordinate, spanned_coverages, 0, "");
                    spanned_blocks.push(new_block);
                    already_processed_coordinates.insert(next_coord.clone());
                }
            }
        }
    }
    spanned_blocks
}

/// Given a coordinate and a list of coordinates in the graph,
/// get coverages for unspanned connections to the prior coordinates.
/// The list of coordinates may include the focus coordinate,
/// but connections between self and self should be omitted.
fn get_unspanned_connections(
    current_coord: &Coordinate,
    previous_coordinates: &[Coordinate],
    clip_coordinates: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
) -> Vec<(ComplexSVBlock, Coordinate)> {
    let mut unspanned_blocks = Vec::new();

    for prev_coord in previous_coordinates.iter() {
        if current_coord == prev_coord {
            continue;
        }
        let mut coord_pair = (prev_coord, current_coord);
        let mut orientation = "+";
        if prev_coord > current_coord {
            coord_pair = (current_coord, prev_coord);
            orientation = "-";
        }
        let prev_coord_opt = Some(coord_pair.0.clone());
        let unspanned_coverages = get_coordinate_alignments_as_coverages(
            clip_coordinates,
            coord_pair.1,
            prev_coord_opt,
            true,
            false,
            false,
        );
        if !unspanned_coverages.is_empty() {
            let spanning_coordinate_opt: Option<Coordinate> =
                coordinate_from_coverages(&unspanned_coverages);
            if let Some(mut unspanned_coordinate) = spanning_coordinate_opt {
                let variant_ids: Vec<String> = coord_pair
                    .1
                    .variant_ids
                    .iter()
                    .chain(prev_coord.variant_ids.iter())
                    .cloned()
                    .collect();
                unspanned_coordinate.variant_ids = variant_ids;

                let new_block =
                    ComplexSVBlock::new(unspanned_coordinate, BTreeMap::new(), 0, orientation);
                unspanned_blocks.push((new_block, current_coord.clone()));
            }
        }
    }

    unspanned_blocks
}

/// Add orientations where possible, from each end toward the middle
fn add_orientations(annotated_event_graph: &mut [ComplexSVBlock]) {
    // if there is a single last coordinate, the end orientation can be inferred as +
    let mut sample_order_idx_counts = HashMap::new();
    let mut max_idx = 0;
    for block in annotated_event_graph.iter() {
        sample_order_idx_counts
            .entry(block.sample_order_index)
            .or_insert(Vec::new())
            .push(block.clone());
        if block.sample_order_index > max_idx {
            max_idx = block.sample_order_index;
        }
    }
    let mut last_known_idx = 0;
    if let Some(final_blocks) = sample_order_idx_counts.get(&max_idx) {
        if final_blocks.len() == 1 {
            if let Some(final_block) = final_blocks.first() {
                if !final_block.coverages.is_empty() {
                    last_known_idx = max_idx;
                }
            }
        }
    }
    // add orientations from the beginning toward the end
    add_orientations_directionally(annotated_event_graph, last_known_idx, false);
    // add orientations from the end toward the beginning
    add_orientations_directionally(annotated_event_graph, last_known_idx, true);
}

/// Identify the orientation of the breakends in the complex event by starting the beginning and
/// inferring the orientation relative to the start in each case where a direct connection exists.
fn add_orientations_directionally(
    annotated_event_graph: &mut [ComplexSVBlock],
    last_known_sample_order_idx: u32,
    reverse_iterate: bool,
) {
    if reverse_iterate {
        annotated_event_graph.reverse();
    }
    let mut prev_block_opt: Option<ComplexSVBlock> = None;
    for block in annotated_event_graph.iter_mut() {
        if block.sample_order_index == 0 || block.sample_order_index == last_known_sample_order_idx
        {
            block.orientation = Orientation::Forward;
            prev_block_opt = Some(block.clone());
        } else if !(block.orientation == Orientation::Missing) {
            // if we've reached a block that already has a known orientation, we're done
            prev_block_opt = Some(block.clone());
        } else if let Some(prev_block) = prev_block_opt {
            if prev_block.orientation == Orientation::Missing {
                break; // once the chain of connections from the beginning has broken, we can't infer orientations anymore
            }
            let inner_orientation_forward = if reverse_iterate {
                Orientation::Reverse
            } else {
                Orientation::Forward
            };
            let inner_orientation_reverse = if reverse_iterate {
                Orientation::Forward
            } else {
                Orientation::Reverse
            };

            if block.coverages.is_empty() && !(prev_block.orientation == Orientation::Missing) {
                // if this block is unspanned and the previous one has a known orientation
                // we can assume this one's order from start and end
                if block.region.start < block.region.end {
                    block.orientation = inner_orientation_forward;
                } else {
                    block.orientation = inner_orientation_reverse;
                }
                prev_block_opt = Some(block.clone());
            } else if prev_block.coverages.is_empty()
                && !(prev_block.orientation == Orientation::Missing)
            {
                // if the previous block has a known orientation and is an
                // unspanned connecting block, we infer the current block's orientation from previous
                let mut prev_block_end =
                    Coordinate::new(prev_block.region.end_chrom.clone(), prev_block.region.end);
                if prev_block.orientation == inner_orientation_reverse {
                    prev_block_end = Coordinate::new(
                        prev_block.region.start_chrom.clone(),
                        prev_block.region.start,
                    );
                }

                let curr_block_start =
                    Coordinate::new(block.region.start_chrom.clone(), block.region.start);
                let curr_block_end =
                    Coordinate::new(block.region.end_chrom.clone(), block.region.end);
                if prev_block_end.is_within(&curr_block_start) {
                    // previous unspanned block connects to the start of the current block,
                    // which means the orientation is forward
                    block.orientation = inner_orientation_forward;
                } else if prev_block_end.is_within(&curr_block_end) {
                    // previous unspanned block connects to the end of the current block,
                    // which means the orientation is reverse
                    block.orientation = inner_orientation_reverse;
                }
                prev_block_opt = Some(block.clone());
            } else {
                // if neither the previous nor current block is an unspanned connecting block,
                // we can't infer any more orientations because we've reached an ambiguous connection
                break;
            }
        } else {
            break;
        }
    }
    if reverse_iterate {
        annotated_event_graph.reverse();
    }
}

/// Create a new coordinate struct to represent the full block spanned by a
/// coverages hashmap
fn coordinate_from_coverages(coverages: &BTreeMap<Coordinate, u32>) -> Option<Coordinate> {
    let leftmost_opt = coverages.keys().min();
    let rightmost_opt = coverages.keys().max();
    if let (Some(leftmost), Some(rightmost)) = (leftmost_opt, rightmost_opt) {
        let new_coord = Coordinate {
            start_chrom: leftmost.start_chrom.clone(),
            start: leftmost.start,
            end_chrom: rightmost.end_chrom.clone(),
            end: rightmost.end,
            confidence_interval: leftmost.confidence_interval,
            variant_ids: Vec::new(),
        };
        Some(new_coord)
    } else {
        None
    }
}

/// If alignments can be found to match a single coordinate on the unclipped side, or to connect a pair of alignments,
/// extract them and convert them to coverage hashtables
fn get_coordinate_alignments_as_coverages(
    clip_alignments: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    first_coordinate: &Coordinate,
    second_coordinate_opt: Option<Coordinate>,
    get_unspanned: bool,
    is_first: bool,
    is_last: bool,
) -> BTreeMap<Coordinate, u32> {
    let mut first_break_alignments = get_first_coordinate_alignments(
        first_coordinate,
        clip_alignments,
        is_first,
        is_last,
        get_unspanned,
    );

    //find alignments that support the second coordinate
    let second_break_alignments;
    if let Some(second_coordinate) = &second_coordinate_opt {
        second_break_alignments = get_second_break_alignments(
            first_coordinate,
            second_coordinate,
            &mut first_break_alignments,
            clip_alignments,
            get_unspanned,
        );
    } else {
        // if there's no second coordinate, just return the coverages from the first
        return alignments_to_coverage(
            &first_break_alignments,
            get_unspanned,
            first_coordinate,
            &None,
        );
    }

    let first_alignments_set: HashSet<FwdStrandSplitReadSegment> =
        first_break_alignments.iter().cloned().collect();
    let second_alignments_set: HashSet<FwdStrandSplitReadSegment> =
        second_break_alignments.into_iter().collect();

    let shared_alignments: Vec<FwdStrandSplitReadSegment> = first_alignments_set
        .intersection(&second_alignments_set)
        .cloned()
        .collect();

    if !shared_alignments.is_empty() {
        let all_alignments: Vec<FwdStrandSplitReadSegment> = first_alignments_set
            .union(&second_alignments_set)
            .cloned()
            .collect();
        // return the coverage from all clipped reads for both coordinates
        alignments_to_coverage(
            &all_alignments,
            get_unspanned,
            first_coordinate,
            &second_coordinate_opt,
        )
    } else {
        // return empty coverage map
        BTreeMap::new()
    }
}

/// Find alignments that match a first coordinate
fn get_first_coordinate_alignments(
    first_coordinate: &Coordinate,
    clip_alignments: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    is_first: bool,
    is_last: bool,
    get_unspanned: bool,
) -> Vec<FwdStrandSplitReadSegment> {
    let mut first_break_alignments = Vec::new();
    let clip_cluster_keys: Vec<&Coordinate> = clip_alignments.keys().collect();
    for clip_cluster in clip_cluster_keys.iter() {
        if first_coordinate.is_within(clip_cluster) {
            if let Some(coord_alignments) = clip_alignments.get(clip_cluster) {
                if get_unspanned {
                    first_break_alignments.append(
                        &mut coord_alignments
                            .iter()
                            .filter(|&a| a.fwd_read_end == a.fwd_read_start)
                            .cloned()
                            .collect(),
                    )
                } else {
                    first_break_alignments.append(
                        &mut coord_alignments
                            .iter()
                            .filter(|a| a.fwd_read_end != a.fwd_read_start)
                            .cloned()
                            .collect(),
                    );
                }

                if is_first {
                    first_break_alignments.retain(|a| {
                        a.end <= (first_coordinate.start + first_coordinate.confidence_interval.0)
                    });
                } else if is_last {
                    first_break_alignments.retain(|a| {
                        a.pos >= (first_coordinate.start - first_coordinate.confidence_interval.0)
                    });
                }
            }
        }
    }
    first_break_alignments
}

/// Find alignments that match a second coordinate
fn get_second_break_alignments(
    first_coordinate: &Coordinate,
    second_coordinate: &Coordinate,
    first_break_alignments: &mut Vec<FwdStrandSplitReadSegment>,
    clip_alignments: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    get_unspanned: bool,
) -> Vec<FwdStrandSplitReadSegment> {
    let mut second_break_alignments = Vec::new();
    let clip_cluster_keys: Vec<&Coordinate> = clip_alignments.keys().collect();
    for clip_cluster in clip_cluster_keys.iter() {
        if second_coordinate.is_within(clip_cluster) {
            let left_coord = std::cmp::min(first_coordinate, &second_coordinate);
            let right_coord = std::cmp::max(first_coordinate, &second_coordinate);
            let strict_left_side = left_coord.start;
            let permissive_left_side = left_coord.start - left_coord.confidence_interval.0;
            let strict_right_side = right_coord.end;
            let permissive_right_side = right_coord.end + right_coord.confidence_interval.0;

            first_break_alignments.retain(|a| {
                ((a.pos >= permissive_left_side)
                    && (a.end <= permissive_right_side)
                    && a.is_start_softclipped
                    && a.is_end_softclipped)
                    || (a.pos >= strict_left_side) && (a.end <= strict_right_side)
            });
            if let Some(coord_alignments) = clip_alignments.get(clip_cluster) {
                let mut second_coord_alignments;
                if get_unspanned {
                    second_coord_alignments = coord_alignments
                        .iter()
                        .filter(|a| !a.spans && is_alignment_in_bounds(a, left_coord, right_coord))
                        .cloned()
                        .collect();
                } else {
                    second_coord_alignments = coord_alignments
                        .iter()
                        .filter(|a| a.spans && is_alignment_in_bounds(a, left_coord, right_coord))
                        .cloned()
                        .collect();
                }
                second_break_alignments.append(&mut second_coord_alignments);
            }
        }
    }
    second_break_alignments
}

/// Function to determine in an alignment is within the bounds of a pair of coordinates
fn is_alignment_in_bounds(
    a: &FwdStrandSplitReadSegment,
    left_coord: &Coordinate,
    right_coord: &Coordinate,
) -> bool {
    let strict_left_side = left_coord.start;
    let permissive_left_side = left_coord.start - left_coord.confidence_interval.0;
    let strict_right_side = right_coord.end;
    let permissive_right_side = right_coord.end + right_coord.confidence_interval.0;

    let permissively_in_bounds =
        (a.pos >= permissive_left_side) && (a.end <= permissive_right_side);
    let double_clipped = a.is_start_softclipped && a.is_end_softclipped;
    let strictly_in_bounds = (a.pos >= strict_left_side) && (a.end <= strict_right_side);

    let alignment_in_bounds = strictly_in_bounds || (permissively_in_bounds && double_clipped);
    alignment_in_bounds
}

/// Use a vector of alignments to generate a BTreeMap of coverages for a block,
/// Where each entry is a coordinate defining a block and the depth of coverage for that block
pub fn alignments_to_coverage(
    alignments: &Vec<FwdStrandSplitReadSegment>,
    get_unspanned: bool,
    start_coord: &Coordinate,
    end_coord_opt: &Option<Coordinate>,
) -> BTreeMap<Coordinate, u32> {
    let mut starts = BTreeMap::new();
    let mut ends = BTreeMap::new();
    let mut breaks = Vec::new();
    let end_coord = match end_coord_opt {
        Some(end_coord) => end_coord,
        None => start_coord,
    };

    for alignment in alignments {
        if (get_unspanned && alignment.spans) || (!get_unspanned && !alignment.spans) {
            continue;
        }
        let start = Coordinate::new(alignment.chrom.clone(), alignment.pos);
        let end = Coordinate::new(alignment.second_chrom.clone(), alignment.end);
        let alignment_start_matches_block_start = (alignment.chrom == start_coord.start_chrom)
            && (alignment.pos.abs_diff(start_coord.start) as i64
                <= start_coord.confidence_interval.0);

        let alignment_start_matches_block_end = (alignment.chrom == end_coord.end_chrom)
            && (alignment.pos.abs_diff(end_coord.start) as i64 <= end_coord.confidence_interval.0);

        let alignment_end_matches_block_start = (alignment.second_chrom == start_coord.start_chrom)
            && (alignment.end.abs_diff(start_coord.start) as i64
                <= start_coord.confidence_interval.0);

        let alignment_end_matches_block_end = (alignment.second_chrom == end_coord.end_chrom)
            && (alignment.end.abs_diff(end_coord.start) as i64 <= end_coord.confidence_interval.0);

        let mut start_clip_matches = (!alignment.is_start_softclipped)
            || (alignment_start_matches_block_start || alignment_start_matches_block_end);
        let mut end_clip_matches = (!alignment.is_end_softclipped)
            || (alignment_end_matches_block_start || alignment_end_matches_block_end);

        // if the start and end match, which will be the case when the
        // alignments are for a single-ended coordinate, they only have to match on one end
        if start_coord == end_coord {
            start_clip_matches = alignment_start_matches_block_start
                || alignment_start_matches_block_end
                || alignment_end_matches_block_start
                || alignment_end_matches_block_end;
            end_clip_matches = start_clip_matches;
        }

        if !start_clip_matches || !end_clip_matches {
            continue;
        }
        if !starts.contains_key(&start) {
            breaks.push(start.clone());
        }
        if !ends.contains_key(&end) {
            breaks.push(end.clone());
        }
        *starts.entry(start).or_insert(0) += 1;
        *ends.entry(end).or_insert(0) += 1;
    }

    breaks.sort();
    let mut coverage = 0;
    let mut coverages = BTreeMap::new();
    for coordinate in breaks {
        if let Some(start_count) = starts.get(&coordinate) {
            coverage += start_count;
        }
        if let Some(end_count) = ends.get(&coordinate) {
            if coverage >= *end_count {
                coverage -= end_count;
            }
        }
        coverages.insert(coordinate, coverage);
    }
    coverages
}
