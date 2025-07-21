use std::{
    collections::{BTreeMap, HashMap, HashSet},
    time::SystemTime,
};

use log::debug;

use crate::{
    containers::{
        AnnotationProcessingState, BlockProcessingContext, BlockProcessingState, ComplexSVBlock,
        ComplexSVCalls, Coordinate, EventGraph, FwdStrandSplitReadSegment, Orientation,
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
                log::debug!("Skipped event: {start} -> {end}",);
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
///    any at the previous level. Add these as unspanned connections.
/// 4) Add spanned: Find coverage for spanned connections going forward from the current coordinate.
/// 5) Add tied connections: Find coverage for unspanned connections between this coordinate
///    and any others at the same sample order index, then add them as unspanned connections.
/// 6) When all levels have been traversed, add a final spanned coordinate from the last one going downstream.
fn annotate_graph(
    event_graph: &HashMap<u32, Vec<Coordinate>>,
    clip_coordinates: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
) -> Vec<ComplexSVBlock> {
    let mut annotated_event_graph: Vec<ComplexSVBlock> = Vec::new();

    if let Some(max_idx) = event_graph.keys().max() {
        let max_idx = *max_idx as i64;
        let mut processing_state = initialize_annotation_state();

        process_all_sample_indices(
            event_graph,
            clip_coordinates,
            &mut annotated_event_graph,
            &mut processing_state,
            max_idx,
        );

        finalize_annotation(
            clip_coordinates,
            &mut annotated_event_graph,
            event_graph,
            max_idx as u32,
        );
    }

    annotated_event_graph
}

/// Initializes the processing state for graph annotation
fn initialize_annotation_state() -> AnnotationProcessingState {
    AnnotationProcessingState {
        already_processed_coordinates: HashSet::new(),
        last_event_spanned: false,
        order_scaler: 0,
    }
}

/// Processes all sample indices in the event graph
fn process_all_sample_indices(
    event_graph: &HashMap<u32, Vec<Coordinate>>,
    clip_coordinates: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    annotated_event_graph: &mut Vec<ComplexSVBlock>,
    processing_state: &mut AnnotationProcessingState,
    max_idx: i64,
) {
    for sample_idx in 0..max_idx + 1 {
        process_single_sample_index(
            event_graph,
            clip_coordinates,
            annotated_event_graph,
            processing_state,
            sample_idx,
            max_idx,
        );
    }
}

/// Processes coordinates for a single sample index
fn process_single_sample_index(
    event_graph: &HashMap<u32, Vec<Coordinate>>,
    clip_coordinates: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    annotated_event_graph: &mut Vec<ComplexSVBlock>,
    processing_state: &mut AnnotationProcessingState,
    sample_idx: i64,
    max_idx: i64,
) {
    if let Some(coordinates) = event_graph.get(&(sample_idx as u32)).cloned() {
        let prev_coords = get_previous_coordinates(event_graph, sample_idx);
        let (mut context, mut state) = setup_processing_context(
            event_graph,
            clip_coordinates,
            annotated_event_graph,
            processing_state,
            sample_idx,
            max_idx,
            prev_coords,
        );

        let (last_event_spanned, order_scaler) = add_blocks(&coordinates, &mut context, &mut state);

        // Update processing state with results
        processing_state.last_event_spanned = last_event_spanned;
        processing_state.order_scaler = order_scaler;
    }
}

/// Gets the previous coordinates for a given sample index
fn get_previous_coordinates(
    event_graph: &HashMap<u32, Vec<Coordinate>>,
    sample_idx: i64,
) -> Vec<Coordinate> {
    let mut prev_coords = Vec::new();
    if sample_idx > 0 {
        if let Some(found_prev_coords) = event_graph.get(&(sample_idx as u32 - 1)) {
            prev_coords.append(&mut found_prev_coords.clone());
        }
    }
    prev_coords
}

/// Sets up the processing context and state for a sample index
fn setup_processing_context<'a>(
    event_graph: &'a HashMap<u32, Vec<Coordinate>>,
    clip_coordinates: &'a HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    annotated_event_graph: &'a mut Vec<ComplexSVBlock>,
    processing_state: &'a mut AnnotationProcessingState,
    sample_idx: i64,
    max_idx: i64,
    prev_coords: Vec<Coordinate>,
) -> (BlockProcessingContext<'a>, BlockProcessingState) {
    let context = BlockProcessingContext {
        event_graph,
        clip_coordinates,
        already_processed_coordinates: &mut processing_state.already_processed_coordinates,
        annotated_event_graph,
    };

    let state = BlockProcessingState {
        sample_idx,
        is_first: sample_idx == 0,
        is_last: sample_idx == max_idx,
        prev_coords,
        last_event_spanned: processing_state.last_event_spanned,
        order_scaler: processing_state.order_scaler,
    };

    (context, state)
}

/// Finalizes the annotation by adding the last block
fn finalize_annotation(
    clip_coordinates: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    annotated_event_graph: &mut Vec<ComplexSVBlock>,
    event_graph: &HashMap<u32, Vec<Coordinate>>,
    max_sample_idx: u32,
) {
    add_last_block(
        clip_coordinates,
        annotated_event_graph,
        event_graph,
        max_sample_idx,
    );
}

/// Loop through and add coordinate blocks that are tied at the same
/// level of the sample order to the annotated event graph.
fn add_blocks(
    coordinates: &[Coordinate],
    context: &mut BlockProcessingContext,
    state: &mut BlockProcessingState,
) -> (bool, u32) {
    let mut max_order_scaler = 0;
    for (i, coordinate) in coordinates.iter().enumerate() {
        let mut curr_order_scaler = state.order_scaler;
        if context.already_processed_coordinates.contains(coordinate) {
            continue;
        }
        if state.is_first {
            add_first_block(
                context.clip_coordinates,
                coordinate,
                context.annotated_event_graph,
                &mut state.prev_coords,
            );
            state.last_event_spanned = true;
        } else {
            let mut unspanned_added = false;
            for (mut unspanned_block, _) in
                get_unspanned_connections(coordinate, &state.prev_coords, context.clip_coordinates)
            {
                unspanned_block.sample_order_index = curr_order_scaler + state.sample_idx as u32;
                context.annotated_event_graph.push(unspanned_block);
                unspanned_added = true;
                context
                    .already_processed_coordinates
                    .insert(coordinate.clone());
                state.last_event_spanned = false;
            }
            let blocks_to_add: Vec<(ComplexSVBlock, Coordinate)> = generate_spanned_blocks(
                state.sample_idx,
                context.event_graph,
                context.clip_coordinates,
                coordinate,
                context.already_processed_coordinates,
                state.is_last,
            );
            if state.last_event_spanned {
                curr_order_scaler += 1;
            }
            for (mut spanned_block, _) in blocks_to_add {
                spanned_block.sample_order_index = curr_order_scaler + state.sample_idx as u32;
                if unspanned_added {
                    spanned_block.sample_order_index += 1;
                }
                context.annotated_event_graph.push(spanned_block);
                state.last_event_spanned = true;
            }
        }
        if curr_order_scaler > max_order_scaler {
            max_order_scaler = curr_order_scaler;
        }

        if coordinates.len() > 1 && i == 0 {
            if let Some(tied_coords) = context.event_graph.get(&(state.sample_idx as u32)) {
                let tied_unspanned_blocks =
                    get_unspanned_connections(coordinate, tied_coords, context.clip_coordinates);
                for (mut tied_unspanned_block, _) in tied_unspanned_blocks {
                    tied_unspanned_block.sample_order_index =
                        curr_order_scaler + state.sample_idx as u32;
                    context.annotated_event_graph.push(tied_unspanned_block);
                }
            };
        }
    }
    (state.last_event_spanned, max_order_scaler)
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
///    perform a fuzzy match check to see if it's been added already (match to previous end coord).
///    If not already added, continue with the following steps.
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

    let orientation_context = OrientationInferenceContext::new(reverse_iterate);
    let mut prev_block_opt: Option<ComplexSVBlock> = None;

    for block in annotated_event_graph.iter_mut() {
        if should_set_forward_orientation(block, last_known_sample_order_idx) {
            set_block_forward_orientation(block, &mut prev_block_opt);
        } else if has_known_orientation(block) {
            prev_block_opt = Some(block.clone());
        } else if let Some(prev_block) = prev_block_opt.clone() {
            if !can_continue_orientation_inference(&prev_block) {
                break;
            }

            if !infer_block_orientation(
                block,
                &prev_block,
                &orientation_context,
                &mut prev_block_opt,
            ) {
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

/// Context for orientation inference based on iteration direction
struct OrientationInferenceContext {
    forward: Orientation,
    reverse: Orientation,
}

impl OrientationInferenceContext {
    fn new(reverse_iterate: bool) -> Self {
        if reverse_iterate {
            Self {
                forward: Orientation::Reverse,
                reverse: Orientation::Forward,
            }
        } else {
            Self {
                forward: Orientation::Forward,
                reverse: Orientation::Reverse,
            }
        }
    }
}

/// Determines if a block should be set to forward orientation
fn should_set_forward_orientation(
    block: &ComplexSVBlock,
    last_known_sample_order_idx: u32,
) -> bool {
    block.sample_order_index == 0 || block.sample_order_index == last_known_sample_order_idx
}

/// Sets a block to forward orientation and updates previous block
fn set_block_forward_orientation(
    block: &mut ComplexSVBlock,
    prev_block_opt: &mut Option<ComplexSVBlock>,
) {
    block.orientation = Orientation::Forward;
    *prev_block_opt = Some(block.clone());
}

/// Checks if a block already has a known orientation
fn has_known_orientation(block: &ComplexSVBlock) -> bool {
    !(block.orientation == Orientation::Missing)
}

/// Determines if orientation inference can continue based on previous block
fn can_continue_orientation_inference(prev_block: &ComplexSVBlock) -> bool {
    prev_block.orientation != Orientation::Missing
}

/// Attempts to infer block orientation based on previous block
fn infer_block_orientation(
    block: &mut ComplexSVBlock,
    prev_block: &ComplexSVBlock,
    context: &OrientationInferenceContext,
    prev_block_opt: &mut Option<ComplexSVBlock>,
) -> bool {
    if is_unspanned_block(block) && has_known_orientation(prev_block) {
        infer_unspanned_block_orientation(block, context, prev_block_opt);
        true
    } else if is_unspanned_block(prev_block) && has_known_orientation(prev_block) {
        infer_from_unspanned_previous_block(block, prev_block, context, prev_block_opt)
    } else {
        // Can't infer more orientations - ambiguous connection
        false
    }
}

/// Checks if a block is unspanned (has no coverages)
fn is_unspanned_block(block: &ComplexSVBlock) -> bool {
    block.coverages.is_empty()
}

/// Infers orientation for an unspanned block based on region coordinates
fn infer_unspanned_block_orientation(
    block: &mut ComplexSVBlock,
    context: &OrientationInferenceContext,
    prev_block_opt: &mut Option<ComplexSVBlock>,
) {
    if block.region.start < block.region.end {
        block.orientation = context.forward;
    } else {
        block.orientation = context.reverse;
    }
    *prev_block_opt = Some(block.clone());
}

/// Infers current block orientation from an unspanned previous block
fn infer_from_unspanned_previous_block(
    block: &mut ComplexSVBlock,
    prev_block: &ComplexSVBlock,
    context: &OrientationInferenceContext,
    prev_block_opt: &mut Option<ComplexSVBlock>,
) -> bool {
    let prev_block_end = get_previous_block_effective_end(prev_block, context);
    let curr_block_start = Coordinate::new(block.region.start_chrom.clone(), block.region.start);
    let curr_block_end = Coordinate::new(block.region.end_chrom.clone(), block.region.end);

    if prev_block_end.is_within(&curr_block_start) {
        block.orientation = context.forward;
        *prev_block_opt = Some(block.clone());
        true
    } else if prev_block_end.is_within(&curr_block_end) {
        block.orientation = context.reverse;
        *prev_block_opt = Some(block.clone());
        true
    } else {
        // Could not determine connection point
        true
    }
}

/// Gets the effective end coordinate of a previous block based on its orientation
fn get_previous_block_effective_end(
    prev_block: &ComplexSVBlock,
    context: &OrientationInferenceContext,
) -> Coordinate {
    if prev_block.orientation == context.reverse {
        Coordinate::new(
            prev_block.region.start_chrom.clone(),
            prev_block.region.start,
        )
    } else {
        Coordinate::new(prev_block.region.end_chrom.clone(), prev_block.region.end)
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
            let left_coord = std::cmp::min(first_coordinate, second_coordinate);
            let right_coord = std::cmp::max(first_coordinate, second_coordinate);
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

    strictly_in_bounds || (permissively_in_bounds && double_clipped)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils;

    fn create_test_coordinate(chrom: &str, pos: i64) -> Coordinate {
        Coordinate::new(chrom.to_string(), pos)
    }

    fn create_test_event_graph() -> HashMap<u32, Vec<Coordinate>> {
        let mut event_graph = HashMap::new();
        event_graph.insert(0, vec![create_test_coordinate("chr1", 1000)]);
        event_graph.insert(1, vec![create_test_coordinate("chr1", 2000)]);
        event_graph.insert(2, vec![create_test_coordinate("chr1", 3000)]);
        event_graph
    }

    fn create_test_clip_coordinates() -> HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>> {
        let mut clip_coordinates = HashMap::new();

        let coord1 = create_test_coordinate("chr1", 1000);
        let coord2 = create_test_coordinate("chr1", 2000);
        let coord3 = create_test_coordinate("chr1", 3000);

        let mut alignments1 = HashSet::new();
        let shared_alignment =
            utils::create_test_alignment_with_clips("chr1", 995, 1010, "read1", true, true);
        alignments1.insert(shared_alignment.clone());
        alignments1.insert(utils::create_test_alignment_with_clips(
            "chr1", 990, 1095, "read2", true, true,
        ));

        let mut alignments2 = HashSet::new();
        alignments2.insert(shared_alignment);
        alignments2.insert(utils::create_test_alignment_with_clips(
            "chr1", 1995, 2100, "read3", true, true,
        ));
        alignments2.insert(utils::create_test_alignment_with_clips(
            "chr1", 1990, 2095, "read4", true, true,
        ));

        let mut alignments3 = HashSet::new();
        alignments3.insert(utils::create_test_alignment_with_clips(
            "chr1", 2995, 3100, "read5", true, true,
        ));
        alignments3.insert(utils::create_test_alignment_with_clips(
            "chr1", 2990, 3095, "read6", true, true,
        ));

        clip_coordinates.insert(coord1, alignments1);
        clip_coordinates.insert(coord2, alignments2);
        clip_coordinates.insert(coord3, alignments3);

        clip_coordinates
    }

    #[test]
    fn test_get_event_ends() {
        let event_graph = create_test_event_graph();
        let (start, end) = get_event_ends(&event_graph);

        assert!(start.is_some());
        assert!(end.is_some());
        assert_eq!(start.unwrap().start, 1000);
        assert_eq!(end.unwrap().start, 3000);
    }

    #[test]
    fn test_get_event_ends_empty() {
        let event_graph = HashMap::new();
        let (start, end) = get_event_ends(&event_graph);

        assert!(start.is_none());
        assert!(end.is_none());
    }

    #[test]
    fn test_initialize_annotation_state() {
        let state = initialize_annotation_state();

        assert!(state.already_processed_coordinates.is_empty());
        assert!(!state.last_event_spanned);
        assert_eq!(state.order_scaler, 0);
    }

    #[test]
    fn test_get_previous_coordinates() {
        let event_graph = create_test_event_graph();

        // Test for index 1 (should return coordinates from index 0)
        let prev_coords = get_previous_coordinates(&event_graph, 1);
        assert_eq!(prev_coords.len(), 1);
        assert_eq!(prev_coords[0].start, 1000);

        // Test for index 0 (should return empty vector)
        let prev_coords = get_previous_coordinates(&event_graph, 0);
        assert!(prev_coords.is_empty());
    }

    #[test]
    fn test_coordinate_from_coverages() {
        let mut coverages = BTreeMap::new();
        coverages.insert(create_test_coordinate("chr1", 1000), 5);
        coverages.insert(create_test_coordinate("chr1", 2000), 3);

        let coord = coordinate_from_coverages(&coverages);
        assert!(coord.is_some());
        assert_eq!(coord.unwrap().start, 1000); // Should return first coordinate
    }

    #[test]
    fn test_coordinate_from_coverages_empty() {
        let coverages = BTreeMap::new();
        let coord = coordinate_from_coverages(&coverages);
        assert!(coord.is_none());
    }

    #[test]
    fn test_is_alignment_in_bounds() {
        let left_coord = create_test_coordinate("chr1", 1000);
        let right_coord = create_test_coordinate("chr1", 2000);

        // Alignment within bounds
        let alignment = utils::create_test_alignment("chr1", 1500, 1600, "read1");
        assert!(is_alignment_in_bounds(
            &alignment,
            &left_coord,
            &right_coord
        ));

        // Alignment outside bounds (before left)
        let alignment = utils::create_test_alignment("chr1", 500, 600, "read2");
        assert!(!is_alignment_in_bounds(
            &alignment,
            &left_coord,
            &right_coord
        ));

        // Alignment outside bounds (after right)
        let alignment = utils::create_test_alignment("chr1", 2500, 2600, "read3");
        assert!(!is_alignment_in_bounds(
            &alignment,
            &left_coord,
            &right_coord
        ));
    }

    #[test]
    fn test_alignments_to_coverage() {
        let start_coord = create_test_coordinate("chr1", 1000);
        let end_coord = Some(create_test_coordinate("chr1", 2000));

        let alignments = vec![
            utils::create_test_alignment("chr1", 995, 1100, "read1"),
            utils::create_test_alignment("chr1", 990, 1095, "read2"),
            utils::create_test_alignment("chr1", 1995, 2100, "read3"),
        ];

        let coverage = alignments_to_coverage(&alignments, false, &start_coord, &end_coord);

        assert!(!coverage.is_empty());
        // Should have coverage at the alignment positions
        assert!(coverage.contains_key(&create_test_coordinate("chr1", 995)));
        assert!(coverage.contains_key(&create_test_coordinate("chr1", 990)));
        assert!(coverage.contains_key(&create_test_coordinate("chr1", 1995)));
    }

    #[test]
    fn test_alignments_to_coverage_empty() {
        let start_coord = create_test_coordinate("chr1", 1000);
        let end_coord = Some(create_test_coordinate("chr1", 2000));
        let alignments = vec![];

        let coverage = alignments_to_coverage(&alignments, false, &start_coord, &end_coord);
        assert!(coverage.is_empty());
    }

    #[test]
    fn test_has_known_orientation() {
        let mut block =
            ComplexSVBlock::new(create_test_coordinate("chr1", 1000), BTreeMap::new(), 0, "");

        // Initially no orientation
        assert!(!has_known_orientation(&block));

        // Set orientation
        block.orientation = Orientation::Forward;
        assert!(has_known_orientation(&block));
    }

    #[test]
    fn test_is_unspanned_block() {
        let mut block = ComplexSVBlock::new(
            create_test_coordinate("chr1", 1000),
            BTreeMap::new(),
            0,
            "+",
        );

        // Initially unspanned (no coverages)
        assert!(is_unspanned_block(&block));

        // Add coverages
        block.coverages.insert("chr1:1000-1100".to_string(), 5);
        assert!(!is_unspanned_block(&block));
    }

    #[test]
    fn test_can_continue_orientation_inference() {
        let mut block =
            ComplexSVBlock::new(create_test_coordinate("chr1", 1000), BTreeMap::new(), 0, "");

        // Block without known orientation
        assert!(!can_continue_orientation_inference(&block));

        // Block with known orientation
        block.orientation = Orientation::Forward;
        assert!(can_continue_orientation_inference(&block));
    }

    #[test]
    fn test_should_set_forward_orientation() {
        let block = ComplexSVBlock::new(
            create_test_coordinate("chr1", 1000),
            BTreeMap::new(),
            5,
            "+",
        );

        // Should set forward orientation when sample order index equals last known
        assert!(should_set_forward_orientation(&block, 5));

        // Should not set forward orientation when sample order index is different
        assert!(!should_set_forward_orientation(&block, 3));
    }

    #[test]
    fn test_annotate_graphs_empty() {
        let event_graphs = vec![];
        let clip_coordinates = HashMap::new();

        let result = annotate_graphs(&event_graphs, &clip_coordinates);
        assert!(result.event_graphs.is_empty());
    }

    #[test]
    fn test_annotate_graphs_single_event() {
        let mut event_graphs = vec![];
        let event_graph = EventGraph {
            graph: create_test_event_graph(),
        };
        event_graphs.push(event_graph);

        let clip_coordinates = create_test_clip_coordinates();

        let result = annotate_graphs(&event_graphs, &clip_coordinates);
        assert!(!result.event_graphs.is_empty());
    }

    #[test]
    fn test_annotate_graphs_large_graph() {
        let mut event_graphs = vec![];
        let mut large_graph = HashMap::new();

        // Create a graph larger than MAX_GRAPH_SIZE
        for i in 0..MAX_GRAPH_SIZE + 10 {
            large_graph.insert(
                i as u32,
                vec![create_test_coordinate("chr1", i as i64 * 1000)],
            );
        }

        let event_graph = EventGraph { graph: large_graph };
        event_graphs.push(event_graph);

        let clip_coordinates = HashMap::new();

        let result = annotate_graphs(&event_graphs, &clip_coordinates);
        // Large graphs should be skipped
        assert!(result.event_graphs.is_empty());
    }

    #[test]
    fn test_annotate_graph() {
        let event_graph = create_test_event_graph();
        let clip_coordinates = create_test_clip_coordinates();

        let result = annotate_graph(&event_graph, &clip_coordinates);
        assert!(!result.is_empty());
    }

    #[test]
    fn test_annotate_graph_empty() {
        let event_graph = HashMap::new();
        let clip_coordinates = HashMap::new();

        let result = annotate_graph(&event_graph, &clip_coordinates);
        assert!(result.is_empty());
    }

    #[test]
    fn test_get_coordinate_alignments_as_coverages() {
        let clip_coordinates = create_test_clip_coordinates();
        let first_coord = create_test_coordinate("chr1", 1000);
        let second_coord = Some(create_test_coordinate("chr1", 2000));

        let coverages = get_coordinate_alignments_as_coverages(
            &clip_coordinates,
            &first_coord,
            second_coord,
            false,
            true,
            false,
        );

        // println!("DEBUG: coverages = {:?}", coverages);
        assert!(!coverages.is_empty());
    }

    #[test]
    fn test_get_first_coordinate_alignments() {
        let clip_coordinates = create_test_clip_coordinates();
        let first_coord = create_test_coordinate("chr1", 1000);

        let alignments =
            get_first_coordinate_alignments(&first_coord, &clip_coordinates, true, false, false);

        assert!(!alignments.is_empty());
    }

    #[test]
    fn test_get_second_break_alignments() {
        let clip_coordinates = create_test_clip_coordinates();
        let first_coord = create_test_coordinate("chr1", 1000);
        let second_coord = create_test_coordinate("chr1", 2000);
        let mut first_break_alignments =
            vec![utils::create_test_alignment("chr1", 995, 1100, "read1")];

        let alignments = get_second_break_alignments(
            &first_coord,
            &second_coord,
            &mut first_break_alignments,
            &clip_coordinates,
            false,
        );

        // Should return alignments that match both coordinates
        assert!(!alignments.is_empty());
    }

    #[test]
    fn test_add_orientations() {
        let mut blocks = vec![
            ComplexSVBlock::new(
                create_test_coordinate("chr1", 1000),
                BTreeMap::new(),
                0,
                "+",
            ),
            ComplexSVBlock::new(
                create_test_coordinate("chr1", 2000),
                BTreeMap::new(),
                1,
                "+",
            ),
            ComplexSVBlock::new(
                create_test_coordinate("chr1", 3000),
                BTreeMap::new(),
                2,
                "+",
            ),
        ];

        // Set some orientations
        blocks[0].orientation = Orientation::Forward;
        blocks[1].sample_order_index = 1;
        blocks[2].sample_order_index = 2;

        add_orientations(&mut blocks);

        // Should have inferred orientations
        assert!(blocks.iter().any(has_known_orientation));
    }

    #[test]
    fn test_orientation_inference_context() {
        let context = OrientationInferenceContext::new(false);
        assert_eq!(context.forward, Orientation::Forward);
        assert_eq!(context.reverse, Orientation::Reverse);

        let context = OrientationInferenceContext::new(true);
        assert_eq!(context.forward, Orientation::Reverse);
        assert_eq!(context.reverse, Orientation::Forward);
    }

    #[test]
    fn test_get_previous_blocks() {
        let mut blocks = vec![
            ComplexSVBlock::new(
                create_test_coordinate("chr1", 1000),
                BTreeMap::new(),
                0,
                "+",
            ),
            ComplexSVBlock::new(
                create_test_coordinate("chr1", 2000),
                BTreeMap::new(),
                1,
                "+",
            ),
            ComplexSVBlock::new(
                create_test_coordinate("chr1", 3000),
                BTreeMap::new(),
                2,
                "+",
            ),
        ];

        let prev_blocks = get_previous_blocks(2, &mut blocks);
        assert!(!prev_blocks.is_empty());
    }

    #[test]
    fn test_already_added() {
        let mut prev_ends = HashSet::new();
        let block = ComplexSVBlock::new(
            create_test_coordinate("chr1", 1000),
            BTreeMap::new(),
            0,
            "+",
        );
        prev_ends.insert(block);

        // Test with coordinate that is NOT within the block's end coordinate
        let coord = create_test_coordinate("chr1", 2000);
        assert!(!already_added(&prev_ends, &coord));

        // Test with coordinate that IS within the block's end coordinate
        let coord_within = create_test_coordinate("chr1", 1005);
        assert!(already_added(&prev_ends, &coord_within));
    }

    #[test]
    fn test_edge_cases() {
        // Test with single coordinate
        let mut event_graph = HashMap::new();
        event_graph.insert(0, vec![create_test_coordinate("chr1", 1000)]);

        let clip_coordinates = create_test_clip_coordinates();
        let result = annotate_graph(&event_graph, &clip_coordinates);

        // Should handle single coordinate case
        assert!(!result.is_empty());
    }

    #[test]
    fn test_coordinate_equality() {
        let coord1 = create_test_coordinate("chr1", 1000);
        let coord2 = create_test_coordinate("chr1", 1000);
        let coord3 = create_test_coordinate("chr1", 2000);

        assert_eq!(coord1, coord2);
        assert_ne!(coord1, coord3);
    }

    #[test]
    fn test_complex_sv_block_creation() {
        let block = ComplexSVBlock::new(
            create_test_coordinate("chr1", 1000),
            BTreeMap::new(),
            0,
            "+",
        );

        assert!(block.coverages.is_empty());
        assert!(has_known_orientation(&block));
        assert_eq!(block.sample_order_index, 0);
    }

    #[test]
    fn test_complex_sv_calls_creation() {
        let blocks = vec![ComplexSVBlock::new(
            create_test_coordinate("chr1", 1000),
            BTreeMap::new(),
            0,
            "+",
        )];
        let event_graphs = vec![blocks];
        let calls = ComplexSVCalls::new(event_graphs);

        assert_eq!(calls.event_graphs.len(), 1);
    }
}
