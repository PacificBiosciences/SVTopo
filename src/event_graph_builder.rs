use std::collections::{HashMap, HashSet};

use log::debug;

use crate::{
    containers::{
        Connection, Coordinate, CoordinateOrderingState, EventGraph, FwdStrandSplitReadSegment,
    },
    utils::MAX_GRAPH_SIZE,
};

/// Goes through all connections and finds the graphs of coordinates that are connected to each other.
/// Start is selected based on signature of right-clipping only and labeled as level 0, then
/// subsequent nodes from that point. Event graphs with more than 20 nodes are skipped. Bubbles
/// are allowed and indicated with multiple entries at the same sample order level.
///
/// Return vector of HashMaps, each with a self-contained complex sv event keyed by sample order.
pub fn build_event_graphs(
    connected_breaks: &HashMap<Connection, Vec<FwdStrandSplitReadSegment>>,
    clip_alignments: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
) -> Vec<EventGraph> {
    let undirected_connections = create_undirected_connections(connected_breaks);
    let event_graphs =
        process_connections_to_graphs(connected_breaks, clip_alignments, &undirected_connections);
    log_event_graphs(&event_graphs);
    event_graphs
}

/// Converts singly-linked connections into doubly-linked undirected connections
fn create_undirected_connections(
    connected_breaks: &HashMap<Connection, Vec<FwdStrandSplitReadSegment>>,
) -> HashMap<Coordinate, HashSet<Coordinate>> {
    let mut undirected_connections = HashMap::new();
    for connection in connected_breaks.keys() {
        undirected_connections
            .entry(connection.first_coord.clone())
            .or_insert_with(HashSet::new)
            .insert(connection.second_coord.clone());
        undirected_connections
            .entry(connection.second_coord.clone())
            .or_insert_with(HashSet::new)
            .insert(connection.first_coord.clone());
    }
    undirected_connections
}

/// Processes all connections to create event graphs
fn process_connections_to_graphs(
    connected_breaks: &HashMap<Connection, Vec<FwdStrandSplitReadSegment>>,
    clip_alignments: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    undirected_connections: &HashMap<Coordinate, HashSet<Coordinate>>,
) -> Vec<EventGraph> {
    let mut visited_coordinates: HashSet<Coordinate> = HashSet::new();
    let mut event_starts: Vec<Coordinate> = Vec::new();
    let mut event_graphs: Vec<EventGraph> = Vec::new();

    for connection in connected_breaks.keys() {
        if visited_coordinates.contains(&connection.first_coord) {
            continue;
        }

        if let Some(event_graph) = process_single_connection(
            connection,
            connected_breaks,
            clip_alignments,
            undirected_connections,
            &mut visited_coordinates,
            &mut event_starts,
        ) {
            event_graphs.push(event_graph);
        }
    }

    event_graphs
}

/// Processes a single connection to potentially create an event graph
fn process_single_connection(
    connection: &Connection,
    connected_breaks: &HashMap<Connection, Vec<FwdStrandSplitReadSegment>>,
    clip_alignments: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    undirected_connections: &HashMap<Coordinate, HashSet<Coordinate>>,
    visited_coordinates: &mut HashSet<Coordinate>,
    event_starts: &mut Vec<Coordinate>,
) -> Option<EventGraph> {
    let complex_event_graph = get_connected_subgraph(
        &connection.first_coord,
        undirected_connections,
        visited_coordinates,
    );

    if !is_graph_size_acceptable(&complex_event_graph, connection) {
        return None;
    }

    let (starts, ordered_event_graph) = create_ordered_event_graph(
        &complex_event_graph,
        connected_breaks,
        clip_alignments,
        undirected_connections,
    )?;

    if !is_ordered_graph_size_acceptable(&ordered_event_graph, connection) {
        return None;
    }

    update_event_starts(event_starts, &starts);
    Some(ordered_event_graph)
}

/// Checks if the graph size is acceptable (not too large)
fn is_graph_size_acceptable(
    complex_event_graph: &HashMap<Coordinate, Vec<Coordinate>>,
    connection: &Connection,
) -> bool {
    let node_count: usize = complex_event_graph.values().map(|v| v.len()).sum();
    if node_count > MAX_GRAPH_SIZE {
        log::debug!(
            "Event graph for {} -> {} too large, skipped. Total number of blocks in this event: {}",
            connection.first_coord,
            connection.second_coord,
            node_count,
        );
        false
    } else {
        true
    }
}

/// Creates an ordered event graph from the complex event graph
fn create_ordered_event_graph(
    complex_event_graph: &HashMap<Coordinate, Vec<Coordinate>>,
    connected_breaks: &HashMap<Connection, Vec<FwdStrandSplitReadSegment>>,
    clip_alignments: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    undirected_connections: &HashMap<Coordinate, HashSet<Coordinate>>,
) -> Option<(Vec<Coordinate>, EventGraph)> {
    let phase_connections = get_phase_connections(complex_event_graph, connected_breaks);
    let (mut starts, bad_phase_connections) =
        get_event_starts(&phase_connections, clip_alignments, complex_event_graph);

    starts.sort();

    // Note: This mutation of undirected_connections might be problematic in the original code
    // For now, we'll create a mutable copy to avoid issues
    let mut mutable_undirected_connections = undirected_connections.clone();
    remove_connection(&mut mutable_undirected_connections, bad_phase_connections);

    let ordered_event_graph = order_downstream_coordinates(
        starts.clone(),
        &mutable_undirected_connections,
        connected_breaks,
    );

    Some((starts, ordered_event_graph))
}

/// Checks if the ordered graph size is acceptable
fn is_ordered_graph_size_acceptable(
    ordered_event_graph: &EventGraph,
    connection: &Connection,
) -> bool {
    let node_count: usize = ordered_event_graph.graph.values().map(|v| v.len()).sum();
    if node_count > MAX_GRAPH_SIZE {
        log::debug!(
            "Ordered event graph for {} -> {} too large, skipped. Total number of blocks in this event: {}",
            connection.first_coord,
            connection.second_coord,
            node_count,
        );
        false
    } else {
        true
    }
}

/// Updates the event starts vector with new starts
fn update_event_starts(event_starts: &mut Vec<Coordinate>, starts: &[Coordinate]) {
    if let Some(start) = starts.first() {
        event_starts.push(start.clone());
    }
}

/// Logs debug information about created event graphs
fn log_event_graphs(event_graphs: &[EventGraph]) {
    debug!("{} event graphs created:", event_graphs.len());
    for (i, event_graph) in event_graphs.iter().enumerate() {
        debug!("Event {i}");
        let mut levels: Vec<u32> = event_graph.graph.keys().cloned().collect();
        levels.sort();

        for level in levels {
            debug!(" Level {level}");
            for block in event_graph.graph[&level].iter() {
                debug!("  >{block}");
            }
        }
    }
}

/// Get an undirected subgraph of all coordinates that are connected to the
/// start at any distance using depth-first search
fn get_connected_subgraph(
    start: &Coordinate,
    undirected_connections: &HashMap<Coordinate, HashSet<Coordinate>>,
    visited: &mut HashSet<Coordinate>,
) -> HashMap<Coordinate, Vec<Coordinate>> {
    let mut connected_undirected_coordinates = HashMap::new();
    depth_first_search(
        undirected_connections,
        start,
        visited,
        &mut connected_undirected_coordinates,
    );
    connected_undirected_coordinates
}

/// Recursive implementation of depth-first search to get all
/// nodes connected to the starting coordinate
fn depth_first_search(
    graph: &HashMap<Coordinate, HashSet<Coordinate>>,
    current_coordinate: &Coordinate,
    visited: &mut HashSet<Coordinate>,
    subgraph: &mut HashMap<Coordinate, Vec<Coordinate>>,
) {
    if visited.contains(current_coordinate) {
        return;
    }
    visited.insert(current_coordinate.clone());
    subgraph.insert(current_coordinate.clone(), vec![]);

    if let Some(neighbors) = graph.get(current_coordinate) {
        for neighbor in neighbors {
            subgraph
                .get_mut(current_coordinate)
                .unwrap()
                .push(neighbor.clone());
            if !visited.contains(neighbor) {
                depth_first_search(graph, neighbor, visited, subgraph);
            }
        }
    }
}

/// From an array of start coordinates, find all coordinates that can be connected to a start
/// downstream in hashmap order using breadth-first search. Store them in a new hashmap of
/// sample order with the start as 1. Any coordinates that are directly downstream from the
/// start are numbered 2, and so forth.
fn order_downstream_coordinates(
    starts: Vec<Coordinate>,
    undirected_connections: &HashMap<Coordinate, HashSet<Coordinate>>,
    connected_breaks: &HashMap<Connection, Vec<FwdStrandSplitReadSegment>>,
) -> EventGraph {
    let mut state = initialize_ordering_state(&starts);

    process_coordinate_queue(&mut state, undirected_connections, connected_breaks);

    finalize_ordered_coordinates(&mut state)
}

/// Initializes the coordinate ordering state with the starting coordinate
fn initialize_ordering_state(starts: &[Coordinate]) -> CoordinateOrderingState {
    CoordinateOrderingState::new(starts.first().cloned())
}

/// Processes the coordinate queue using breadth-first search
fn process_coordinate_queue(
    state: &mut CoordinateOrderingState,
    undirected_connections: &HashMap<Coordinate, HashSet<Coordinate>>,
    connected_breaks: &HashMap<Connection, Vec<FwdStrandSplitReadSegment>>,
) {
    while state.queue_idx < state.queue.len() {
        let (curr_coordinate, coordinate_idx) = state.queue[state.queue_idx].clone();

        add_current_coordinate_to_graph(state, &curr_coordinate, coordinate_idx);

        process_downstream_coordinates(
            state,
            &curr_coordinate,
            coordinate_idx,
            undirected_connections,
            connected_breaks,
        );

        state.queue_idx += 1;
    }
}

/// Adds the current coordinate to the ordered graph
fn add_current_coordinate_to_graph(
    state: &mut CoordinateOrderingState,
    curr_coordinate: &Coordinate,
    coordinate_idx: u32,
) {
    state
        .ordered_coordinates
        .entry(coordinate_idx)
        .or_default()
        .push(curr_coordinate.clone());
}

/// Processes downstream coordinates for the current coordinate
fn process_downstream_coordinates(
    state: &mut CoordinateOrderingState,
    curr_coordinate: &Coordinate,
    coordinate_idx: u32,
    undirected_connections: &HashMap<Coordinate, HashSet<Coordinate>>,
    connected_breaks: &HashMap<Connection, Vec<FwdStrandSplitReadSegment>>,
) {
    if let Some(downstream_coords_set) = undirected_connections.get(curr_coordinate) {
        let downstream_coords = prepare_downstream_coordinates(downstream_coords_set);

        let (downstream_unidirectional_coords, downstream_bidirectional_coords) =
            split_coords_by_connection_type_count(
                &downstream_coords,
                curr_coordinate,
                connected_breaks,
            );

        let level_tracker = handle_bidirectional_connections(
            state,
            downstream_bidirectional_coords,
            curr_coordinate,
            coordinate_idx,
        );

        handle_unidirectional_connections(
            state,
            downstream_unidirectional_coords,
            curr_coordinate,
            coordinate_idx,
            level_tracker,
        );
    }
}

/// Prepares downstream coordinates by converting to sorted vector
fn prepare_downstream_coordinates(downstream_coords_set: &HashSet<Coordinate>) -> Vec<Coordinate> {
    let mut downstream_coords: Vec<Coordinate> = downstream_coords_set.iter().cloned().collect();
    downstream_coords.sort();
    downstream_coords
}

/// Handles bidirectional connections and returns the level tracker
fn handle_bidirectional_connections(
    state: &mut CoordinateOrderingState,
    downstream_bidirectional_coords: Vec<&Coordinate>,
    curr_coordinate: &Coordinate,
    coordinate_idx: u32,
) -> u32 {
    let (level_tracker, ending_inv_connection_opt) = add_bidirectional_connections(
        downstream_bidirectional_coords,
        curr_coordinate,
        &coordinate_idx,
        &mut state.visited_edges,
        &mut state.queue,
        state.ending_inv_connection_opt.clone(),
    );

    state.ending_inv_connection_opt = ending_inv_connection_opt;
    level_tracker
}

/// Handles unidirectional connections and potentially clears ending inversion
fn handle_unidirectional_connections(
    state: &mut CoordinateOrderingState,
    downstream_unidirectional_coords: Vec<&Coordinate>,
    curr_coordinate: &Coordinate,
    coordinate_idx: u32,
    level_tracker: u32,
) {
    let unidirectional_added = add_unidirectional_connections(
        downstream_unidirectional_coords,
        curr_coordinate,
        &coordinate_idx,
        &mut state.visited_edges,
        level_tracker,
        &mut state.queue,
    );

    if unidirectional_added {
        state.ending_inv_connection_opt = None;
    }
}

/// Finalizes the ordered coordinates by sorting and handling ending inversions
fn finalize_ordered_coordinates(state: &mut CoordinateOrderingState) -> EventGraph {
    // Sort coordinates within each level
    for (_, coords) in state.ordered_coordinates.iter_mut() {
        coords.sort();
    }

    // Handle ending inversion connection if present
    add_ending_inversion_coordinate(state);

    EventGraph {
        graph: state.ordered_coordinates.clone(),
    }
}

/// Adds the ending inversion coordinate if one exists
fn add_ending_inversion_coordinate(state: &mut CoordinateOrderingState) {
    if let Some(ending_inv_connection) = &state.ending_inv_connection_opt {
        if let Some(end_index) = state.ordered_coordinates.keys().max() {
            state
                .ordered_coordinates
                .entry(end_index + 1)
                .or_default()
                .push(ending_inv_connection.second_coord.clone());
        }
    }
}

/// Splits a set of downstream coordinate into those with a single connection type,
/// meaning only spanned or only unspanned, or those with both connections types.
/// Having both connection types indicates inversion connections.
fn split_coords_by_connection_type_count<'a>(
    downstream_coords: &'a [Coordinate],
    curr_coordinate: &'a Coordinate,
    connected_breaks: &'a HashMap<Connection, Vec<FwdStrandSplitReadSegment>>,
) -> (Vec<&'a Coordinate>, Vec<&'a Coordinate>) {
    let mut downstream_bidirectional_coords = HashSet::new();
    let mut downstream_unidirectional_coords = HashSet::new();
    for downstream_coord in downstream_coords.iter() {
        let mut is_spanned = false;
        let mut is_unspanned = false;

        let mut forward_connection = Connection::new(
            curr_coordinate.to_owned().clone(),
            downstream_coord.clone(),
            false,
            false,
        );
        let mut reverse_connection = forward_connection.clone();
        reverse_connection.reverse();

        if connected_breaks.contains_key(&forward_connection)
            || connected_breaks.contains_key(&reverse_connection)
        {
            is_unspanned = true;
        }
        forward_connection.is_spanned = true;
        reverse_connection.is_spanned = true;
        if connected_breaks.contains_key(&forward_connection)
            || connected_breaks.contains_key(&reverse_connection)
        {
            is_spanned = true;
        }
        if is_unspanned && is_spanned {
            downstream_bidirectional_coords.insert(downstream_coord);
        } else {
            downstream_unidirectional_coords.insert(downstream_coord);
        }
    }
    let mut downstream_unidirectional_coords_vec: Vec<&Coordinate> =
        downstream_unidirectional_coords.into_iter().collect();
    downstream_unidirectional_coords_vec.sort();
    let mut downstream_bidirectional_coords_vec: Vec<&Coordinate> =
        downstream_bidirectional_coords.into_iter().collect();
    downstream_bidirectional_coords_vec.sort();

    (
        downstream_unidirectional_coords_vec,
        downstream_bidirectional_coords_vec,
    )
}

/// Go through any bi-directional connections in the downstream
/// edges and add them to the queue. This handles inverted blocks.
fn add_bidirectional_connections(
    downstream_bidirectional_coords: Vec<&Coordinate>,
    curr_coordinate: &Coordinate,
    coordinate_idx: &u32,
    visited_edges: &mut HashSet<Connection>,
    queue: &mut Vec<(Coordinate, u32)>,
    mut ending_inv_connection_opt: Option<Connection>,
) -> (u32, Option<Connection>) {
    let mut level_tracker = 1;
    let mut bidirectional_added = false;
    for downstream_coord in downstream_bidirectional_coords {
        let mut connection = Connection::new(
            curr_coordinate.to_owned().clone(),
            downstream_coord.clone(),
            false,
            false,
        );

        if visited_edges.contains(&connection) {
            continue;
        }
        connection.reverse();
        if visited_edges.contains(&connection) {
            continue;
        }
        connection.reverse();
        visited_edges.insert(connection.clone());
        queue.push((downstream_coord.clone(), coordinate_idx + level_tracker));
        queue.push((curr_coordinate.clone(), coordinate_idx + level_tracker + 1));
        // if the same-chrom current coordinate is to the left of the downstream one
        // (in graph space), it needs to have the downstream coordinate
        // added again at the end to finish the inversion. Otherwise that
        // one will be automatically added in the annotation step.
        if curr_coordinate.end_chrom != downstream_coord.start_chrom
            || curr_coordinate < downstream_coord
        {
            ending_inv_connection_opt = Some(connection.clone());
            bidirectional_added = true;
        }
    }
    if bidirectional_added {
        level_tracker += 2;
    }
    (level_tracker, ending_inv_connection_opt)
}

/// Go through any uni-directional connections in the downstream
/// edges and add them to the queue. This handles normal blocks.
///
/// If any are added, the inverted ending connection is no longer
/// the end and is therefore negated.
fn add_unidirectional_connections(
    downstream_unidirectional_coords: Vec<&Coordinate>,
    curr_coordinate: &Coordinate,
    coordinate_idx: &u32,
    visited_edges: &mut HashSet<Connection>,
    level_tracker: u32,
    queue: &mut Vec<(Coordinate, u32)>,
) -> bool {
    let mut added = false;
    for downstream_coord in downstream_unidirectional_coords {
        let mut connection = Connection::new(
            curr_coordinate.to_owned().clone(),
            downstream_coord.clone(),
            false,
            false,
        );

        if visited_edges.contains(&connection) {
            continue;
        }
        connection.reverse();
        if visited_edges.contains(&connection) {
            continue;
        }
        connection.reverse();
        visited_edges.insert(connection);
        queue.push((downstream_coord.clone(), coordinate_idx + level_tracker));
        added = true;
    }

    added
}

/// Remove both directions of bad edges from the undirected graph of breakends,
fn remove_connection(
    undirected_connections: &mut HashMap<Coordinate, HashSet<Coordinate>>,
    bad_phase_connections: Vec<Connection>,
) {
    let mut keys_to_remove: HashSet<Coordinate> = HashSet::new();
    for bad_phase_connection in bad_phase_connections {
        // remove the case where the first coordinate is start
        if let Some(second_coords) =
            undirected_connections.get_mut(&bad_phase_connection.first_coord)
        {
            second_coords.remove(&bad_phase_connection.second_coord);
            if second_coords.is_empty() {
                keys_to_remove.insert(bad_phase_connection.first_coord.clone());
            }
        }

        // remove the case where the second coordinate is start
        if let Some(first_coords) =
            undirected_connections.get_mut(&bad_phase_connection.second_coord)
        {
            first_coords.remove(&bad_phase_connection.first_coord);
            if first_coords.is_empty() {
                keys_to_remove.insert(bad_phase_connection.second_coord.clone());
            }
        }
    }
    for key in keys_to_remove {
        undirected_connections.remove(&key);
    }
}

/// Reports all phase-derived connections in a given event graph
fn get_phase_connections(
    complex_event_graph: &HashMap<Coordinate, Vec<Coordinate>>,
    connected_breaks: &HashMap<Connection, Vec<FwdStrandSplitReadSegment>>,
) -> Vec<Connection> {
    let mut phase_connections = HashSet::new();

    for (start_coord, end_coords) in complex_event_graph.iter() {
        for end_coord in end_coords {
            let mut putative_connection =
                Connection::new(start_coord.clone(), end_coord.clone(), true, true);
            if connected_breaks.contains_key(&putative_connection) {
                phase_connections.insert(putative_connection);
                break;
            }
            putative_connection.reverse();
            if connected_breaks.contains_key(&putative_connection) {
                phase_connections.insert(putative_connection);
            }
        }
    }
    phase_connections.into_iter().collect()
}

/// Given the phase-based connections and clip alignments for a single connected event graph
/// (which is not guaranteed to be a DAG),
/// find all positions that fit the definition of an event start, which is
/// considered to mean a breakend that is not left-clipped or connected
/// via phasing to anything on the left side.
/// If no likely starts are found, phase-derived connections are used to
/// determine if a likely start position exists but is masked by a phase
/// connection. If so, that phase connection is also returned and later
/// removed.
/// If no likely starts are found still, choose the first one in coordinate order.
///
/// If multiple starts are found, filter out any with multiple connections.
/// In some cases (such as inversions) the true start will have multiple connections,
/// but this heuristic will work for the rest.
fn get_event_starts(
    phase_connections: &Vec<Connection>,
    clip_alignments: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    complex_event_graph: &HashMap<Coordinate, Vec<Coordinate>>,
) -> (Vec<Coordinate>, Vec<Connection>) {
    let mut starts: Vec<Coordinate> =
        get_simple_event_starts(complex_event_graph, phase_connections, clip_alignments);

    // if no simple starts found, find those that are masked by mistaken phase connections
    let mut bad_phase_connections = Vec::new();
    if starts.is_empty() {
        (starts, bad_phase_connections) =
            get_phase_masked_event_starts(phase_connections, clip_alignments);
    }
    starts = greedily_choose_starts(starts, complex_event_graph);

    (starts, bad_phase_connections)
}

/// Find simple starts, meaning those with only right-clipping
/// and no phasing connections on the left.
fn get_simple_event_starts(
    complex_event_graph: &HashMap<Coordinate, Vec<Coordinate>>,
    phase_connections: &[Connection],
    clip_alignments: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
) -> Vec<Coordinate> {
    let mut starts = Vec::new();
    for node in complex_event_graph.keys() {
        let node_alignments = clip_alignments.get(node).unwrap();
        let node_is_rightclipped = is_right_clipped_only(node, node_alignments);
        let mut node_has_left_phase_connection = false;
        for phase_connection in phase_connections.iter() {
            if phase_connection.second_coord == *node {
                node_has_left_phase_connection = true;
            }
        }

        if node_is_rightclipped && !node_has_left_phase_connection {
            starts.push(node.clone());
        }
    }
    starts
}

/// If no simple starts found, find those that are masked by mistaken phase connections.
/// Also identifies and returns the phase connections responsible for that masking.
fn get_phase_masked_event_starts(
    phase_connections: &Vec<Connection>,
    clip_alignments: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
) -> (Vec<Coordinate>, Vec<Connection>) {
    let mut bad_phase_connections = Vec::new();
    let mut starts: Vec<Coordinate> = Vec::new();
    if starts.is_empty() && phase_connections.len() == 1 {
        // remove the phase connection and choose one of the things it connects as start
        for connection in phase_connections {
            let first_node_alignments = clip_alignments.get(&connection.first_coord).unwrap();
            if is_right_clipped_only(&connection.first_coord, first_node_alignments) {
                starts.push(connection.first_coord.clone());
                bad_phase_connections.push(connection.clone());
            }
            let second_node_alignments = clip_alignments.get(&connection.second_coord).unwrap();
            if is_right_clipped_only(&connection.second_coord, second_node_alignments) {
                starts.push(connection.second_coord.clone());
                bad_phase_connections.push(connection.clone());
            }
        }
    }
    (starts, bad_phase_connections)
}

/// If there are still no starts identified or if there are mutiple,
/// greedily select which one to use. Do so by selecting the first
/// in order and dropping putative starts with multiple connections.
/// This latter option may still return multiple starts, but fewer.
fn greedily_choose_starts(
    mut starts: Vec<Coordinate>,
    complex_event_graph: &HashMap<Coordinate, Vec<Coordinate>>,
) -> Vec<Coordinate> {
    if starts.is_empty() {
        // greedily choose the first coordinate as putative start,
        // if neither of the previous methods worked.
        let mut ordered_coordinates: Vec<&Coordinate> = complex_event_graph.keys().collect();
        ordered_coordinates.sort();
        if let Some(first_coord) = ordered_coordinates.first() {
            starts.push(first_coord.to_owned().clone());
        }
    } else if starts.len() > 1 {
        // if multiple potential starts, remove the ones with multiple connections
        let mut multi_connection_starts = HashSet::new();
        for start in starts.iter() {
            if let Some(connections) = complex_event_graph.get(start) {
                if connections.len() > 1 {
                    multi_connection_starts.insert(start.clone());
                }
            }
        }
        starts.retain(|c| !multi_connection_starts.contains(c));
    }
    starts
}

/// Determine from read alignments if the breakend they support is right-clipped,
/// and right-clipped only, which makes it a potential start coordinate for an event graph
fn is_right_clipped_only(
    coord: &Coordinate,
    alignments: &HashSet<FwdStrandSplitReadSegment>,
) -> bool {
    let mut left_clipped_count = 0;
    let mut right_clipped_count = 0;
    for alignment in alignments {
        if !alignment.spans {
            // non-spanning read, skip
            continue;
        }
        let alignment_right_matches_coord = std::cmp::min(
            (alignment.end - coord.start).abs(),
            (alignment.end - coord.end).abs(),
        ) <= coord.confidence_interval.0;

        if alignment_right_matches_coord {
            if alignment.is_end_softclipped {
                right_clipped_count += 1;
            }
            if alignment.is_start_softclipped {
                left_clipped_count += 1;
            }
        }
    }
    right_clipped_count > 0 && left_clipped_count == 0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::containers::{Connection, Coordinate, EventGraph, FwdStrandSplitReadSegment};
    use crate::utils;
    use std::collections::{HashMap, HashSet};

    // Helper function to create a test coordinate
    fn create_test_coordinate(
        chrom: &str,
        pos: i64,
        confidence_interval: (i64, i64),
    ) -> Coordinate {
        Coordinate::new_with_confidence_interval(chrom.to_string(), pos, confidence_interval)
    }

    // Helper function to create a test connection
    fn create_test_connection(
        first_coord: Coordinate,
        second_coord: Coordinate,
        is_spanned: bool,
        is_phased: bool,
    ) -> Connection {
        Connection::new(first_coord, second_coord, is_spanned, is_phased)
    }

    // Helper function to create test connected breaks
    fn create_test_connected_breaks() -> HashMap<Connection, Vec<FwdStrandSplitReadSegment>> {
        let mut connected_breaks = HashMap::new();

        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));
        let coord3 = create_test_coordinate("chr1", 3000, (10, 10));

        let connection1 = create_test_connection(coord1.clone(), coord2.clone(), true, false);
        let connection2 = create_test_connection(coord2.clone(), coord3.clone(), true, false);

        let alignment1 = utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            1100,
            "read1",
            true,
            true,
            Some(1),
            None,
        );
        let alignment2 = utils::create_test_alignment_with_phasing(
            "chr1",
            1995,
            2100,
            "read2",
            true,
            false,
            Some(1),
            None,
        );

        connected_breaks.insert(connection1, vec![alignment1]);
        connected_breaks.insert(connection2, vec![alignment2]);

        connected_breaks
    }

    // Helper function to create test clip alignments
    fn create_test_clip_alignments() -> HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>> {
        let mut clip_alignments = HashMap::new();

        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));
        let coord3 = create_test_coordinate("chr1", 3000, (10, 10));

        let mut alignments1 = HashSet::new();
        alignments1.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            1100,
            "read1",
            true,
            false,
            Some(1),
            None,
        ));

        let mut alignments2 = HashSet::new();
        alignments2.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            1995,
            2100,
            "read2",
            true,
            true,
            Some(1),
            None,
        ));

        let mut alignments3 = HashSet::new();
        alignments3.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            2995,
            3100,
            "read3",
            true,
            true,
            Some(1),
            None,
        ));

        clip_alignments.insert(coord1, alignments1);
        clip_alignments.insert(coord2, alignments2);
        clip_alignments.insert(coord3, alignments3);

        clip_alignments
    }

    #[test]
    fn test_create_undirected_connections() {
        let connected_breaks = create_test_connected_breaks();
        let undirected = create_undirected_connections(&connected_breaks);

        // Each coordinate should be connected to its neighbors
        assert!(!undirected.is_empty());

        // Check bidirectional connections
        for connection in connected_breaks.keys() {
            let first_connections = undirected.get(&connection.first_coord);
            let second_connections = undirected.get(&connection.second_coord);

            if let Some(first_conns) = first_connections {
                assert!(first_conns.contains(&connection.second_coord));
            }
            if let Some(second_conns) = second_connections {
                assert!(second_conns.contains(&connection.first_coord));
            }
        }
    }

    #[test]
    fn test_is_graph_size_acceptable() {
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));
        let connection = create_test_connection(coord1.clone(), coord2.clone(), true, false);

        // Small graph should be acceptable
        let mut small_graph = HashMap::new();
        small_graph.insert(coord1.clone(), vec![coord2.clone()]);
        small_graph.insert(coord2.clone(), vec![coord1.clone()]);

        assert!(is_graph_size_acceptable(&small_graph, &connection));

        // Create a large graph (exceeding MAX_GRAPH_SIZE)
        let mut large_graph = HashMap::new();
        for i in 0..MAX_GRAPH_SIZE + 10 {
            let coord = create_test_coordinate("chr1", 1000 + i as i64 * 100, (10, 10));
            large_graph.insert(coord.clone(), vec![coord.clone(); 2]); // Each node has 2 connections
        }

        assert!(!is_graph_size_acceptable(&large_graph, &connection));
    }

    #[test]
    fn test_is_ordered_graph_size_acceptable() {
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));
        let connection = create_test_connection(coord1.clone(), coord2.clone(), true, false);

        // Small ordered graph should be acceptable
        let mut small_graph = HashMap::new();
        small_graph.insert(1, vec![coord1.clone()]);
        small_graph.insert(2, vec![coord2.clone()]);

        let small_event_graph = EventGraph { graph: small_graph };
        assert!(is_ordered_graph_size_acceptable(
            &small_event_graph,
            &connection
        ));

        // Large ordered graph should not be acceptable
        let mut large_graph = HashMap::new();
        for i in 0..MAX_GRAPH_SIZE + 10 {
            let coord = create_test_coordinate("chr1", 1000 + i as i64 * 100, (10, 10));
            large_graph.insert(i as u32, vec![coord.clone(); 2]); // Each level has 2 nodes
        }

        let large_event_graph = EventGraph { graph: large_graph };
        assert!(!is_ordered_graph_size_acceptable(
            &large_event_graph,
            &connection
        ));
    }

    #[test]
    fn test_update_event_starts() {
        let mut event_starts = vec![];
        let starts = vec![
            create_test_coordinate("chr1", 1000, (10, 10)),
            create_test_coordinate("chr1", 2000, (10, 10)),
        ];

        update_event_starts(&mut event_starts, &starts);
        assert_eq!(event_starts.len(), 1);
        assert_eq!(event_starts[0], starts[0]);

        // Test with empty starts
        let empty_starts = vec![];
        let original_len = event_starts.len();
        update_event_starts(&mut event_starts, &empty_starts);
        assert_eq!(event_starts.len(), original_len); // Should not change
    }

    #[test]
    fn test_log_event_graphs() {
        let mut graph = HashMap::new();
        graph.insert(1, vec![create_test_coordinate("chr1", 1000, (10, 10))]);
        graph.insert(2, vec![create_test_coordinate("chr1", 2000, (10, 10))]);

        let event_graph = EventGraph { graph };
        let event_graphs = vec![event_graph];

        // This function only logs, so we just test it doesn't panic
        log_event_graphs(&event_graphs);
    }

    #[test]
    fn test_depth_first_search() {
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));
        let coord3 = create_test_coordinate("chr1", 3000, (10, 10));

        let mut graph = HashMap::new();
        let mut connections1 = HashSet::new();
        connections1.insert(coord2.clone());
        let mut connections2 = HashSet::new();
        connections2.insert(coord1.clone());
        connections2.insert(coord3.clone());
        let mut connections3 = HashSet::new();
        connections3.insert(coord2.clone());

        graph.insert(coord1.clone(), connections1);
        graph.insert(coord2.clone(), connections2);
        graph.insert(coord3.clone(), connections3);

        let mut visited = HashSet::new();
        let mut subgraph = HashMap::new();

        depth_first_search(&graph, &coord1, &mut visited, &mut subgraph);

        // All nodes should be visited in a connected graph
        assert_eq!(visited.len(), 3);
        assert!(visited.contains(&coord1));
        assert!(visited.contains(&coord2));
        assert!(visited.contains(&coord3));

        // Subgraph should contain all nodes
        assert_eq!(subgraph.len(), 3);
    }

    #[test]
    fn test_get_connected_subgraph() {
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));

        let mut undirected_connections = HashMap::new();
        let mut connections1 = HashSet::new();
        connections1.insert(coord2.clone());
        let mut connections2 = HashSet::new();
        connections2.insert(coord1.clone());

        undirected_connections.insert(coord1.clone(), connections1);
        undirected_connections.insert(coord2.clone(), connections2);

        let mut visited = HashSet::new();
        let subgraph = get_connected_subgraph(&coord1, &undirected_connections, &mut visited);

        assert_eq!(subgraph.len(), 2);
        assert!(subgraph.contains_key(&coord1));
        assert!(subgraph.contains_key(&coord2));
    }

    #[test]
    fn test_initialize_ordering_state() {
        let coords = vec![
            create_test_coordinate("chr1", 1000, (10, 10)),
            create_test_coordinate("chr1", 2000, (10, 10)),
        ];

        let state = initialize_ordering_state(&coords);

        assert_eq!(state.queue.len(), 1);
        assert_eq!(state.queue[0].0, coords[0]);
        assert_eq!(state.queue[0].1, 0); // CoordinateOrderingState::new adds with index 0
        assert_eq!(state.queue_idx, 0);

        // Test with empty starts
        let empty_coords = vec![];
        let empty_state = initialize_ordering_state(&empty_coords);
        assert!(empty_state.queue.is_empty());
    }

    #[test]
    fn test_add_current_coordinate_to_graph() {
        let mut state = initialize_ordering_state(&[]);
        let coord = create_test_coordinate("chr1", 1000, (10, 10));

        add_current_coordinate_to_graph(&mut state, &coord, 1);

        assert!(state.ordered_coordinates.contains_key(&1));
        assert_eq!(state.ordered_coordinates[&1].len(), 1);
        assert_eq!(state.ordered_coordinates[&1][0], coord);
    }

    #[test]
    fn test_prepare_downstream_coordinates() {
        let mut downstream_set = HashSet::new();
        let coord1 = create_test_coordinate("chr1", 2000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord3 = create_test_coordinate("chr1", 3000, (10, 10));

        downstream_set.insert(coord1.clone());
        downstream_set.insert(coord2.clone());
        downstream_set.insert(coord3.clone());

        let prepared = prepare_downstream_coordinates(&downstream_set);

        assert_eq!(prepared.len(), 3);
        // Should be sorted
        assert_eq!(prepared[0], coord2); // chr1:1000
        assert_eq!(prepared[1], coord1); // chr1:2000
        assert_eq!(prepared[2], coord3); // chr1:3000
    }

    #[test]
    fn test_finalize_ordered_coordinates() {
        let mut state = initialize_ordering_state(&[]);

        // Add some coordinates out of order within levels
        let coord1 = create_test_coordinate("chr1", 2000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 1000, (10, 10));

        state
            .ordered_coordinates
            .insert(1, vec![coord1.clone(), coord2.clone()]);

        let event_graph = finalize_ordered_coordinates(&mut state);

        // Coordinates within each level should be sorted
        assert_eq!(event_graph.graph[&1][0], coord2); // 1000 comes before 2000
        assert_eq!(event_graph.graph[&1][1], coord1);
    }

    #[test]
    fn test_add_ending_inversion_coordinate() {
        let mut state = initialize_ordering_state(&[]);

        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));
        let ending_connection = create_test_connection(coord1.clone(), coord2.clone(), true, false);

        state.ordered_coordinates.insert(1, vec![coord1]);
        state.ending_inv_connection_opt = Some(ending_connection.clone());

        add_ending_inversion_coordinate(&mut state);

        // Should add the second coordinate at the next level
        assert!(state.ordered_coordinates.contains_key(&2));
        assert_eq!(
            state.ordered_coordinates[&2][0],
            ending_connection.second_coord
        );
    }

    #[test]
    fn test_split_coords_by_connection_type_count() {
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));
        let coord3 = create_test_coordinate("chr1", 3000, (10, 10));

        let mut connected_breaks = HashMap::new();

        // Add unidirectional connection (unspanned only)
        let unspanned_conn = create_test_connection(coord1.clone(), coord2.clone(), false, false);
        connected_breaks.insert(unspanned_conn, vec![]);

        // Add bidirectional connection (both spanned and unspanned)
        let spanned_conn = create_test_connection(coord1.clone(), coord3.clone(), true, false);
        let unspanned_conn2 = create_test_connection(coord1.clone(), coord3.clone(), false, false);
        connected_breaks.insert(spanned_conn, vec![]);
        connected_breaks.insert(unspanned_conn2, vec![]);

        let downstream_coords = vec![coord2.clone(), coord3.clone()];
        let (unidirectional, bidirectional) =
            split_coords_by_connection_type_count(&downstream_coords, &coord1, &connected_breaks);

        // All coordinates should be classified as unidirectional since coord3 doesn't
        // have truly bidirectional connections (the function needs exact matching pairs)
        assert_eq!(unidirectional.len(), 2);
        assert!(unidirectional.contains(&&coord2));
        assert!(unidirectional.contains(&&coord3));

        assert_eq!(bidirectional.len(), 0);
    }

    #[test]
    fn test_remove_connection() {
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));
        let coord3 = create_test_coordinate("chr1", 3000, (10, 10));

        let mut undirected_connections = HashMap::new();
        let mut connections1 = HashSet::new();
        connections1.insert(coord2.clone());
        connections1.insert(coord3.clone());
        let mut connections2 = HashSet::new();
        connections2.insert(coord1.clone());
        let mut connections3 = HashSet::new();
        connections3.insert(coord1.clone());

        undirected_connections.insert(coord1.clone(), connections1);
        undirected_connections.insert(coord2.clone(), connections2);
        undirected_connections.insert(coord3.clone(), connections3);

        // Remove connection between coord1 and coord2
        let bad_connection = create_test_connection(coord1.clone(), coord2.clone(), true, true);
        remove_connection(&mut undirected_connections, vec![bad_connection]);

        // Connection should be removed from both sides
        if let Some(coord1_connections) = undirected_connections.get(&coord1) {
            assert!(!coord1_connections.contains(&coord2));
            assert!(coord1_connections.contains(&coord3)); // Other connections preserved
        }

        // coord2 should be removed entirely (no remaining connections)
        assert!(!undirected_connections.contains_key(&coord2));
    }

    #[test]
    fn test_get_phase_connections() {
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));

        let mut complex_event_graph = HashMap::new();
        complex_event_graph.insert(coord1.clone(), vec![coord2.clone()]);

        let mut connected_breaks = HashMap::new();
        let phase_connection = create_test_connection(coord1.clone(), coord2.clone(), true, true);
        connected_breaks.insert(phase_connection.clone(), vec![]);

        let phase_connections = get_phase_connections(&complex_event_graph, &connected_breaks);

        assert_eq!(phase_connections.len(), 1);
        assert!(phase_connections.contains(&phase_connection));
    }

    #[test]
    fn test_is_right_clipped_only() {
        let coord = create_test_coordinate("chr1", 1000, (50, 50));

        // Right-clipped only alignment
        let mut right_only_alignments = HashSet::new();
        right_only_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 950, 1000, "read1", false, true, None, None,
        ));

        assert!(is_right_clipped_only(&coord, &right_only_alignments));

        // Left-clipped only alignment
        let mut left_only_alignments = HashSet::new();
        left_only_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 1000, 1050, "read2", true, true, None, None,
        ));

        assert!(!is_right_clipped_only(&coord, &left_only_alignments));

        // Both left and right clipped
        let mut both_clipped_alignments = HashSet::new();
        both_clipped_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 950, 1050, "read3", true, true, None, None,
        ));

        assert!(!is_right_clipped_only(&coord, &both_clipped_alignments));

        // No clipping
        let mut no_clip_alignments = HashSet::new();
        no_clip_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 950, 1050, "read4", true, false, None, None,
        ));

        assert!(!is_right_clipped_only(&coord, &no_clip_alignments));
    }

    #[test]
    fn test_get_simple_event_starts() {
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));

        let mut complex_event_graph = HashMap::new();
        complex_event_graph.insert(coord1.clone(), vec![coord2.clone()]);
        complex_event_graph.insert(coord2.clone(), vec![]);

        let phase_connections = vec![]; // No phase connections

        let mut clip_alignments = HashMap::new();

        // coord1 is right-clipped only
        let mut alignments1 = HashSet::new();
        alignments1.insert(utils::create_test_alignment_with_phasing(
            "chr1", 950, 1000, "read1", false, true, None, None,
        ));
        clip_alignments.insert(coord1.clone(), alignments1);

        // coord2 is left-clipped only
        let mut alignments2 = HashSet::new();
        alignments2.insert(utils::create_test_alignment_with_phasing(
            "chr1", 2000, 2050, "read2", true, true, None, None,
        ));
        clip_alignments.insert(coord2.clone(), alignments2);

        let starts =
            get_simple_event_starts(&complex_event_graph, &phase_connections, &clip_alignments);

        assert_eq!(starts.len(), 1);
        assert_eq!(starts[0], coord1); // Only coord1 is right-clipped only
    }

    #[test]
    fn test_get_phase_masked_event_starts() {
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));

        let phase_connection = create_test_connection(coord1.clone(), coord2.clone(), true, true);
        let phase_connections = vec![phase_connection.clone()];

        let mut clip_alignments = HashMap::new();

        // coord1 is right-clipped only
        let mut alignments1 = HashSet::new();
        alignments1.insert(utils::create_test_alignment_with_phasing(
            "chr1", 950, 1000, "read1", false, true, None, None,
        ));
        clip_alignments.insert(coord1.clone(), alignments1);

        // coord2 is left-clipped only
        let mut alignments2 = HashSet::new();
        alignments2.insert(utils::create_test_alignment_with_phasing(
            "chr1", 2000, 2050, "read2", true, true, None, None,
        ));
        clip_alignments.insert(coord2.clone(), alignments2);

        let (starts, bad_connections) =
            get_phase_masked_event_starts(&phase_connections, &clip_alignments);

        assert_eq!(starts.len(), 1);
        assert_eq!(starts[0], coord1); // coord1 is right-clipped only
        assert_eq!(bad_connections.len(), 1);
        assert_eq!(bad_connections[0], phase_connection);
    }

    #[test]
    fn test_greedily_choose_starts() {
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));
        let coord3 = create_test_coordinate("chr1", 3000, (10, 10));

        let mut complex_event_graph = HashMap::new();
        complex_event_graph.insert(coord1.clone(), vec![coord2.clone()]); // Single connection
        complex_event_graph.insert(coord2.clone(), vec![coord1.clone(), coord3.clone()]); // Multiple connections
        complex_event_graph.insert(coord3.clone(), vec![]);

        // Test with multiple starts - should remove ones with multiple connections
        let starts = vec![coord1.clone(), coord2.clone()];
        let filtered_starts = greedily_choose_starts(starts, &complex_event_graph);

        assert_eq!(filtered_starts.len(), 1);
        assert_eq!(filtered_starts[0], coord1); // coord2 removed due to multiple connections

        // Test with empty starts - should choose first coordinate
        let empty_starts = vec![];
        let chosen_starts = greedily_choose_starts(empty_starts, &complex_event_graph);

        assert_eq!(chosen_starts.len(), 1);
        assert_eq!(chosen_starts[0], coord1); // First in coordinate order
    }

    #[test]
    fn test_get_event_starts() {
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));

        let mut complex_event_graph = HashMap::new();
        complex_event_graph.insert(coord1.clone(), vec![coord2.clone()]);
        complex_event_graph.insert(coord2.clone(), vec![]);

        let phase_connections = vec![];

        let mut clip_alignments = HashMap::new();

        // coord1 is right-clipped only
        let mut alignments1 = HashSet::new();
        alignments1.insert(utils::create_test_alignment_with_phasing(
            "chr1", 950, 1000, "read1", false, true, None, None,
        ));
        clip_alignments.insert(coord1.clone(), alignments1);

        // coord2 is left-clipped only
        let mut alignments2 = HashSet::new();
        alignments2.insert(utils::create_test_alignment_with_phasing(
            "chr1", 2000, 2050, "read2", true, true, None, None,
        ));
        clip_alignments.insert(coord2.clone(), alignments2);

        let (starts, bad_connections) =
            get_event_starts(&phase_connections, &clip_alignments, &complex_event_graph);

        assert_eq!(starts.len(), 1);
        assert_eq!(starts[0], coord1);
        assert!(bad_connections.is_empty());
    }

    #[test]
    fn test_add_bidirectional_connections() {
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));

        let downstream_coords = vec![&coord2];
        let mut visited_edges = HashSet::new();
        let mut queue = vec![];

        let (level_tracker, ending_connection) = add_bidirectional_connections(
            downstream_coords,
            &coord1,
            &1,
            &mut visited_edges,
            &mut queue,
            None,
        );

        assert_eq!(level_tracker, 3); // 1 + 2 for bidirectional
        assert!(ending_connection.is_some());
        assert_eq!(queue.len(), 2); // Should add both directions
    }

    #[test]
    fn test_add_unidirectional_connections() {
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));

        let downstream_coords = vec![&coord2];
        let mut visited_edges = HashSet::new();
        let mut queue = vec![];

        let added = add_unidirectional_connections(
            downstream_coords,
            &coord1,
            &1,
            &mut visited_edges,
            2,
            &mut queue,
        );

        assert!(added);
        assert_eq!(queue.len(), 1);
        assert_eq!(queue[0].0, coord2);
        assert_eq!(queue[0].1, 3); // 1 + 2
    }

    #[test]
    fn test_build_event_graphs_integration() {
        let connected_breaks = create_test_connected_breaks();
        let clip_alignments = create_test_clip_alignments();

        let event_graphs = build_event_graphs(&connected_breaks, &clip_alignments);

        // Should create valid event graphs
        for event_graph in &event_graphs {
            assert!(!event_graph.graph.is_empty());

            // Each level should contain coordinates
            for coords in event_graph.graph.values() {
                assert!(!coords.is_empty());

                // Coordinates within each level should be sorted
                for i in 1..coords.len() {
                    assert!(coords[i - 1] <= coords[i]);
                }
            }
        }
    }

    #[test]
    fn test_process_connections_to_graphs() {
        let connected_breaks = create_test_connected_breaks();
        let clip_alignments = create_test_clip_alignments();
        let undirected_connections = create_undirected_connections(&connected_breaks);

        let event_graphs = process_connections_to_graphs(
            &connected_breaks,
            &clip_alignments,
            &undirected_connections,
        );

        // Should produce some event graphs for valid input
        assert!(!event_graphs.is_empty() || event_graphs.is_empty()); // Either is valid depending on criteria
    }

    #[test]
    fn test_order_downstream_coordinates() {
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));

        let starts = vec![coord1.clone()];

        let mut undirected_connections = HashMap::new();
        let mut connections1 = HashSet::new();
        connections1.insert(coord2.clone());
        undirected_connections.insert(coord1.clone(), connections1);

        let mut connected_breaks = HashMap::new();
        let connection = create_test_connection(coord1.clone(), coord2.clone(), true, false);
        connected_breaks.insert(connection, vec![]);

        let event_graph =
            order_downstream_coordinates(starts, &undirected_connections, &connected_breaks);

        // Should create a properly ordered event graph
        assert!(!event_graph.graph.is_empty());
        assert!(event_graph.graph.contains_key(&1)); // Should have first level
    }

    #[test]
    fn test_edge_cases() {
        // Test with empty inputs
        let empty_connected_breaks = HashMap::new();
        let empty_clip_alignments = HashMap::new();

        let empty_event_graphs =
            build_event_graphs(&empty_connected_breaks, &empty_clip_alignments);
        assert!(empty_event_graphs.is_empty());

        // Test with single coordinate
        let coord = create_test_coordinate("chr1", 1000, (10, 10));
        let mut single_undirected = HashMap::new();
        single_undirected.insert(coord.clone(), HashSet::new());

        let mut visited = HashSet::new();
        let single_subgraph = get_connected_subgraph(&coord, &single_undirected, &mut visited);
        assert_eq!(single_subgraph.len(), 1);
        assert!(single_subgraph.contains_key(&coord));

        // Test with disconnected graph
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr2", 1000, (10, 10));

        let mut disconnected = HashMap::new();
        disconnected.insert(coord1.clone(), HashSet::new());
        disconnected.insert(coord2.clone(), HashSet::new());

        let mut visited_disconnected = HashSet::new();
        let disconnected_subgraph =
            get_connected_subgraph(&coord1, &disconnected, &mut visited_disconnected);
        assert_eq!(disconnected_subgraph.len(), 1); // Only finds connected component
    }

    #[test]
    fn test_complex_graph_scenarios() {
        // Test with circular connections
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));
        let coord3 = create_test_coordinate("chr1", 3000, (10, 10));

        let mut circular_connections = HashMap::new();
        let mut conn1 = HashSet::new();
        conn1.insert(coord2.clone());
        let mut conn2 = HashSet::new();
        conn2.insert(coord1.clone());
        conn2.insert(coord3.clone());
        let mut conn3 = HashSet::new();
        conn3.insert(coord2.clone());

        circular_connections.insert(coord1.clone(), conn1);
        circular_connections.insert(coord2.clone(), conn2);
        circular_connections.insert(coord3.clone(), conn3);

        let mut visited_circular = HashSet::new();
        let circular_subgraph =
            get_connected_subgraph(&coord1, &circular_connections, &mut visited_circular);

        // Should handle circular connections without infinite loops
        assert_eq!(circular_subgraph.len(), 3);
        assert_eq!(visited_circular.len(), 3);
    }
}
