use std::collections::{HashMap, HashSet};

use log::debug;

use crate::{
    containers::{Connection, Coordinate, EventGraph, FwdStrandSplitReadSegment},
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
    // turn the singly-linked connections into doubly-linked
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

    // Loop through all connections. For each that has not already been incorporated into a complex event graph,
    // traverse the hashmap to extract the connected events in both directions. If a start node can be identified
    // in this subgraph, the entire graph is then ordered and stored.
    let mut visited_coordinates: HashSet<Coordinate> = HashSet::new();
    let mut event_starts: Vec<Coordinate> = Vec::new();
    let mut event_graphs: Vec<EventGraph> = Vec::new();
    for connection in connected_breaks.keys() {
        if visited_coordinates.contains(&connection.first_coord) {
            continue;
        }
        let complex_event_graph = get_connected_subgraph(
            &connection.first_coord,
            &undirected_connections,
            &mut visited_coordinates,
        );

        let node_count: usize = complex_event_graph.values().map(|v| v.len()).sum();
        if node_count > MAX_GRAPH_SIZE {
            log::debug!(
                "Event graph for {} -> {} too large, skipped. Total number of blocks in this event: {}",
                connection.first_coord,
                connection.second_coord,
                node_count,
            );
            continue;
        }

        let phase_connections = get_phase_connections(&complex_event_graph, connected_breaks);

        let (mut starts, bad_phase_connections) =
            get_event_starts(&phase_connections, clip_alignments, &complex_event_graph);
        starts.sort();
        if let Some(start) = starts.first() {
            event_starts.push(start.clone());
        }
        remove_connection(&mut undirected_connections, bad_phase_connections);

        // tracks coordinate order in this complex event
        let ordered_event_graph =
            order_downstream_coordinates(starts, &undirected_connections, connected_breaks);
        let node_count: usize = ordered_event_graph.graph.values().map(|v| v.len()).sum();
        if node_count > MAX_GRAPH_SIZE {
            log::debug!(
                "Ordered event graph for {} -> {} too large, skipped. Total number of blocks in this event: {}",
                connection.first_coord,
                connection.second_coord,
                node_count,
            );
            continue;
        }

        event_graphs.push(ordered_event_graph);
    }

    // section to log output for debug
    debug!("{} event graphs created:", event_graphs.len());
    for (i, event_graph) in event_graphs.iter().enumerate() {
        debug!("Event {}", i);
        let mut levels: Vec<u32> = event_graph.graph.keys().cloned().collect();
        levels.sort();

        for level in levels {
            debug!(" Level {}", level);
            for block in event_graph.graph[&level].iter() {
                debug!("  >{}", block);
            }
        }
    }
    event_graphs
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
    let mut visited_edges = HashSet::new();
    let mut ordered_coordinates = EventGraph::new();
    let mut queue = Vec::new();
    let mut queue_idx = 0;

    if let Some(start) = starts.first() {
        queue.push((start.clone(), 0));
    }

    let mut ending_inv_connection_opt = None;
    while queue_idx < queue.len() {
        let (curr_coordinate, coordinate_idx) = &queue[queue_idx].clone();

        // add the current coordinate to the graph
        ordered_coordinates
            .graph
            .entry(*coordinate_idx)
            .or_default()
            .push((*curr_coordinate).clone());

        // add neighbors of the current coordinate to the queue
        if let Some(downstream_coords_set) = undirected_connections.get(curr_coordinate) {
            let mut downstream_coords: Vec<Coordinate> =
                downstream_coords_set.iter().cloned().collect();
            downstream_coords.sort();

            let (downstream_unidirectional_coords, downstream_bidirectional_coords) =
                split_coords_by_connection_type_count(
                    &downstream_coords,
                    curr_coordinate,
                    connected_breaks,
                );
            let level_tracker;
            (level_tracker, ending_inv_connection_opt) = add_bidirectional_connections(
                downstream_bidirectional_coords,
                curr_coordinate,
                coordinate_idx,
                &mut visited_edges,
                &mut queue,
                ending_inv_connection_opt,
            );

            if add_unidirectional_connections(
                downstream_unidirectional_coords,
                curr_coordinate,
                coordinate_idx,
                &mut visited_edges,
                level_tracker,
                &mut queue,
            ) {
                ending_inv_connection_opt = None;
            }
        }
        queue_idx += 1;
    }

    for (_, coords) in ordered_coordinates.graph.iter_mut() {
        coords.sort();
    }

    if let Some(ending_inv_connection) = ending_inv_connection_opt {
        if let Some(end_index) = ordered_coordinates.graph.keys().max() {
            // add the final INV coordinate to the graph
            ordered_coordinates
                .graph
                .entry(end_index + 1)
                .or_default()
                .push((ending_inv_connection.second_coord.clone()).clone());
        }
    }
    ordered_coordinates
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
