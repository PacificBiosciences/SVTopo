use log::debug;
use std::collections::{HashMap, HashSet};
use std::time::SystemTime;

use crate::containers::{Connection, Coordinate, FwdStrandSplitReadSegment};
use crate::utils;

/// Connect clusters of genomic breaks using alignments, vcf connections, and then phasing
///
/// Returns a hashmap of two connected coordinates and a bool flag indicating alignment support,
/// keyed to a vector of supporting alignments
pub fn connect_clusters(
    cluster_alignments: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    vcf_connections: &Vec<Connection>,
    allow_unphased: bool,
) -> HashMap<Connection, Vec<FwdStrandSplitReadSegment>> {
    let start_time = SystemTime::now();
    let mut non_unique_connections =
        connect_clusters_by_alignment(cluster_alignments, allow_unphased);
    connect_clusters_by_vcf_connection(
        cluster_alignments,
        &mut non_unique_connections,
        vcf_connections,
    );
    connect_clusters_by_phaseset(cluster_alignments, &mut non_unique_connections);

    let mut connections: HashMap<Connection, Vec<FwdStrandSplitReadSegment>> = HashMap::new();
    for (connection, alignments) in non_unique_connections {
        let spanning_tags: HashSet<bool> = alignments.iter().map(|a| a.spans).collect();
        for spanning_tag in spanning_tags {
            let reverse_connection = Connection::new(
                connection.second_coord.clone(),
                connection.first_coord.clone(),
                connection.inferred_from_phasing,
                spanning_tag,
            );
            if let Some(unique_connection_alignments) = connections.get_mut(&reverse_connection) {
                let connection_alignments: Vec<FwdStrandSplitReadSegment> = alignments
                    .iter()
                    .filter(|a| a.spans == spanning_tag)
                    .cloned()
                    .collect();
                unique_connection_alignments.extend(connection_alignments);
                let mut seen = HashSet::new();
                unique_connection_alignments.retain(|a| seen.insert(a.clone()));
            } else {
                connections.insert(connection.clone(), alignments.clone());
            }
        }
    }

    let break_count = connections.len();
    if break_count == 1 {
        debug!("{} connected pair of breaks found", break_count);
    } else {
        debug!("{} connected pairs of breaks found", break_count);
    }

    let mut keys: Vec<&Connection> = connections.keys().collect();
    keys.sort();
    for con in keys {
        debug!(
            "{} -> {} (spans: {}, phased: {})",
            con.first_coord, con.second_coord, con.is_spanned, con.inferred_from_phasing
        );
    }

    debug!(
        "Cluster connection: {}s",
        start_time.elapsed().unwrap().as_secs()
    );
    connections
}

/// Connect clusters of genomic breaks using alignments that are shared between them.
///
/// Returns a hashmap of two connected coordinates and a bool flag indicating alignment support,
/// keyed to a vector of supporting alignments
pub fn connect_clusters_by_alignment(
    cluster_alignments: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    allow_unphased: bool,
) -> HashMap<Connection, Vec<FwdStrandSplitReadSegment>> {
    let mut connections: HashMap<Connection, Vec<FwdStrandSplitReadSegment>> = HashMap::new();

    // find and store all alignments that connect clip clusters
    for (first_cluster, first_cluster_alignments) in cluster_alignments.iter() {
        for (second_cluster, second_cluster_alignments) in cluster_alignments.iter() {
            if first_cluster == second_cluster {
                continue;
            }

            let shared_alignments: HashSet<FwdStrandSplitReadSegment> = first_cluster_alignments
                .intersection(second_cluster_alignments)
                .cloned()
                .collect();
            let spanned_shared_alignments: Vec<FwdStrandSplitReadSegment> = shared_alignments
                .iter()
                .filter(|a| a.spans)
                .cloned()
                .collect();
            let unspanned_shared_alignments: Vec<FwdStrandSplitReadSegment> = shared_alignments
                .iter()
                .filter(|a| !a.spans)
                .cloned()
                .collect();
            let mut connection =
                Connection::new(first_cluster.clone(), second_cluster.clone(), false, true);
            if first_cluster > second_cluster {
                connection.reverse();
            }
            if spanned_shared_alignments.len() >= utils::MIN_CONNECTION_SUPPORT as usize
                && phase_is_consistent(&spanned_shared_alignments, allow_unphased)
            {
                connection.is_spanned = true;
                connections.insert(connection.clone(), spanned_shared_alignments);
            }

            if unspanned_shared_alignments.len() >= utils::MIN_CONNECTION_SUPPORT as usize {
                connection.is_spanned = false;
                connections.insert(connection, unspanned_shared_alignments);
            }
        }
    }

    debug!("{} alignment connections added", connections.len());
    connections
}

/// Use vcf connections where available to add breakend connections.
///
/// Add vcf connections to the hashmap of connections to alignments, with all shared alignments
/// between the two coordinates in the vcf connections included. Spanned and unspanned
/// alignments are stored independently to capture cases where the same two breaks are
/// connected in both ways (as in inversions).
pub fn connect_clusters_by_vcf_connection(
    cluster_alignments: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    connected_breaks: &mut HashMap<Connection, Vec<FwdStrandSplitReadSegment>>,
    vcf_connections: &Vec<Connection>,
) {
    for vcf_connection in vcf_connections {
        if let Some(first_coord_alignments) = cluster_alignments.get(&vcf_connection.first_coord) {
            if let Some(second_coord_alignments) =
                cluster_alignments.get(&vcf_connection.second_coord)
            {
                let shared_alignments: HashSet<FwdStrandSplitReadSegment> = first_coord_alignments
                    .intersection(second_coord_alignments)
                    .cloned()
                    .collect();
                let mut spanned_shared_alignments: HashSet<FwdStrandSplitReadSegment> =
                    HashSet::new();
                let mut unspanned_shared_alignments: HashSet<FwdStrandSplitReadSegment> =
                    HashSet::new();
                for alignment in shared_alignments.iter() {
                    if alignment.spans {
                        spanned_shared_alignments.insert(alignment.clone());
                    } else {
                        unspanned_shared_alignments.insert(alignment.clone());
                    }
                }
                for connection_alignments_set in
                    [spanned_shared_alignments, unspanned_shared_alignments]
                {
                    if connection_alignments_set.len() >= (utils::MIN_CONNECTION_SUPPORT as usize) {
                        let alignments: Vec<FwdStrandSplitReadSegment> =
                            connection_alignments_set.into_iter().collect();
                        let mut connection = vcf_connection.clone();
                        connection.is_spanned = alignments[0].spans;
                        if connection.first_coord > connection.second_coord {
                            connection.reverse();
                        }
                        connected_breaks.insert(connection, alignments);
                    }
                }
            }
        }
    }
}

/// Use phaseset information where available to add breakend connections
///
/// If the left of a breakend is clipped and the right is not, and the reads supporting it are phased,
/// attempt to extend downstream from the break to the next break. If a break can be found where
/// the right is clipped and the left is not, and supporting reads have the phase phaseset ID,
/// assume it's the same block and add a connection. As soon as a new phaseset is found, the chromosome
/// ends, or a megabase of separation is reached, stop looking for an end break to match the start break.
pub fn connect_clusters_by_phaseset(
    cluster_alignments: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    connected_breaks: &mut HashMap<Connection, Vec<FwdStrandSplitReadSegment>>,
) {
    let mut phased_connection_count = 0;
    // a chromosome-keyed map of sorted coordinates
    let mut ordered_breaks: HashMap<String, Vec<&Coordinate>> = HashMap::new();
    for coordinate in cluster_alignments.keys() {
        ordered_breaks
            .entry(coordinate.start_chrom.clone())
            .or_default()
            .push(coordinate);
    }
    for coordinates in ordered_breaks.values_mut() {
        coordinates.sort_unstable();
        for (first_idx, &first_coord) in coordinates.iter().enumerate() {
            let first_coord_alignments = cluster_alignments.get(first_coord).unwrap();
            if is_potential_block_start(first_coord, first_coord_alignments) {
                let mut first_phaseset = 0;
                let mut first_haplotype = 0;
                let mut first_coord_aln_start = i64::MAX;
                let mut first_coord_aln_end = 0;
                for aln in first_coord_alignments {
                    // skip unspanned
                    if !aln.spans {
                        continue;
                    }
                    if aln.pos < first_coord_aln_start {
                        first_coord_aln_start = aln.pos;
                    }
                    if aln.end > first_coord_aln_end {
                        first_coord_aln_end = aln.end;
                    }
                    if let (Some(phaseset), Some(haplotype)) =
                        (&aln.phaseset_tag, &aln.haplotype_tag)
                    {
                        first_phaseset = *phaseset;
                        first_haplotype = *haplotype;
                    }
                }
                for second_coord in coordinates.iter().skip(first_idx + 1) {
                    // assume not a match if not the same chromosome or more than a MAX_PHASE_DIST apart
                    let dist = (first_coord.start - second_coord.start).abs();
                    if second_coord.start_chrom != first_coord.start_chrom
                        || dist > utils::MAX_PHASE_DIST
                    {
                        continue;
                    }
                    // if the alignment support for the first position overlaps the breakpoint for the
                    // second one, the two should have explicit support rather than just phasing support
                    if second_coord.start >= first_coord_aln_start
                        && second_coord.start <= first_coord_aln_end
                    {
                        continue;
                    }
                    let second_coord_alignments = cluster_alignments.get(second_coord).unwrap();
                    let mut second_is_phased = false;
                    let mut second_phase_matches = false;
                    // if any alignments are phased but none match the phase for the first coord,
                    // we've moved to a new phaseset and won't find a match
                    for aln in second_coord_alignments {
                        if let (Some(second_phaseset), Some(second_haplotype)) =
                            (&aln.phaseset_tag, &aln.haplotype_tag)
                        {
                            second_is_phased = true;
                            if *second_phaseset != first_phaseset {
                                break; //we've moved on to a new phaseset
                            } else if *second_haplotype == first_haplotype {
                                second_phase_matches = true;
                                break;
                            }
                        }
                    }
                    if second_is_phased && !second_phase_matches {
                        continue;
                    }
                    if is_block_end(
                        second_coord,
                        second_coord_alignments,
                        first_phaseset,
                        first_haplotype,
                    ) {
                        let combined_alignments = get_combined_alignments_for_phased_connections(
                            first_coord,
                            second_coord,
                            first_coord_alignments,
                            second_coord_alignments,
                        );
                        if combined_alignments.is_empty() {
                            continue;
                        }
                        let spanning_tags: HashSet<bool> =
                            combined_alignments.iter().map(|a| a.spans).collect();
                        for spanning_tag in spanning_tags {
                            let mut new_connection = Connection::new(
                                first_coord.clone(),
                                (*second_coord).clone(),
                                true,
                                spanning_tag,
                            );
                            let connection_alignments: Vec<FwdStrandSplitReadSegment> =
                                combined_alignments
                                    .iter()
                                    .filter(|a| a.spans == spanning_tag)
                                    .cloned()
                                    .collect();
                            if new_connection.first_coord > new_connection.second_coord {
                                new_connection.reverse();
                            }
                            connected_breaks.insert(new_connection, connection_alignments);
                            phased_connection_count += 1;
                        }
                        break;
                    }
                }
            }
        }
    }
    debug!("{} phased connections added", phased_connection_count);
}

/// Find shared alignments that belong to both partners in a phased connection
/// Drop any that are softclipped in a location that doesn't match the coordinates,
/// including a max_clust_dist amount of allowed error.
fn get_combined_alignments_for_phased_connections(
    first_coord: &Coordinate,
    second_coord: &Coordinate,
    first_coord_alignments: &HashSet<FwdStrandSplitReadSegment>,
    second_coord_alignments: &HashSet<FwdStrandSplitReadSegment>,
) -> Vec<FwdStrandSplitReadSegment> {
    let mut shared_alignments: Vec<FwdStrandSplitReadSegment> = Vec::new();
    let bnd_breaks: Vec<(String, i64, i64)> = vec![
        (
            first_coord.start_chrom.clone(),
            first_coord.start,
            first_coord.confidence_interval.0,
        ),
        (
            first_coord.end_chrom.clone(),
            first_coord.end,
            first_coord.confidence_interval.1,
        ),
        (
            second_coord.start_chrom.clone(),
            second_coord.start,
            second_coord.confidence_interval.0,
        ),
        (
            second_coord.end_chrom.clone(),
            second_coord.end,
            second_coord.confidence_interval.1,
        ),
    ];
    for alignment_set in &[first_coord_alignments, second_coord_alignments] {
        for alignment in *alignment_set {
            let mut alignment_start_matches = false;
            let mut alignment_end_matches = false;

            for (break_chrom, break_start, break_ci) in bnd_breaks.iter() {
                if *break_chrom == alignment.chrom {
                    if (break_start - alignment.pos).abs() < *break_ci {
                        alignment_start_matches = true;
                    }
                    if (break_start - alignment.end).abs() < *break_ci {
                        alignment_end_matches = true;
                    }
                }
            }

            if (alignment.is_end_softclipped && !alignment_end_matches)
                || (alignment.is_start_softclipped && !alignment_start_matches)
            {
                continue;
            } else {
                shared_alignments.push(alignment.clone());
            }
        }
    }
    shared_alignments
}

/// Checks alignments that support a coordinate and determines if they are consistently
/// phased to the same phase set, left-clipped, and not right-clipped. If so, the coordinate is
/// a potential start to a block
fn is_potential_block_start(
    coord: &Coordinate,
    alignments: &HashSet<FwdStrandSplitReadSegment>,
) -> bool {
    let mut is_phased = false;
    let mut is_left_clipped = false;
    let mut is_right_clipped = false;

    let mut phase_sets = HashSet::new();
    let mut left_clip_count = 0;
    let mut right_clip_count = 0;
    for alignment in alignments {
        if !alignment.spans || !alignment_matches_break(alignment, coord) {
            continue; // ignore unspanned alignments
        }
        if let Some(phaseset) = &alignment.phaseset_tag {
            phase_sets.insert(phaseset);
        }
        if alignment.is_start_softclipped {
            left_clip_count += 1;
        }
        if alignment.is_end_softclipped {
            right_clip_count += 1;
        }
    }
    if phase_sets.len() == 1 {
        is_phased = true;
    }
    if left_clip_count >= 2 {
        is_left_clipped = true;
    }
    if right_clip_count >= 2 {
        is_right_clipped = true;
    }

    is_phased && is_left_clipped && !is_right_clipped
}

/// Checks alignments that support a coordinate and determines if they are consistently
/// phased to the same phase set as the one passed in,
/// right-clipped, and not left-clipped. If so, the coordinate is
/// the end to a block
fn is_block_end(
    coord: &Coordinate,
    alignments: &HashSet<FwdStrandSplitReadSegment>,
    target_phaseset: i32,
    target_haplotype: i32,
) -> bool {
    let mut phase_matches = false;
    let mut is_left_clipped = false;
    let mut is_right_clipped = false;

    let mut phase_sets = HashSet::new();
    let mut haplotypes = HashSet::new();
    let mut left_clip_count = 0;
    let mut right_clip_count = 0;
    for alignment in alignments {
        if !alignment.spans || !alignment_matches_break(alignment, coord) {
            continue; // ignore unspanned alignments
        }

        if let Some(phaseset) = &alignment.phaseset_tag {
            phase_sets.insert(*phaseset);
        }
        if let Some(haplotype) = &alignment.haplotype_tag {
            haplotypes.insert(*haplotype);
        }
        if alignment.is_start_softclipped {
            left_clip_count += 1;
        }
        if alignment.is_end_softclipped {
            right_clip_count += 1;
        }
    }
    if phase_sets.len() == 1 && haplotypes.len() == 1 {
        if let (Some(aln_phaseset), Some(aln_haplotype)) =
            (phase_sets.iter().next(), haplotypes.iter().next())
        {
            if *aln_phaseset == target_phaseset && *aln_haplotype == target_haplotype {
                phase_matches = true;
            }
        }
    }
    if left_clip_count >= 2 {
        is_left_clipped = true;
    }
    if right_clip_count >= 2 {
        is_right_clipped = true;
    }

    phase_matches && !is_left_clipped && is_right_clipped
}

/// Check if the location of the alignment clip site is
/// within max_clust_dist of the coordinate
fn alignment_matches_break(alignment: &FwdStrandSplitReadSegment, coord: &Coordinate) -> bool {
    let mut matches = false;
    if (coord.start - coord.confidence_interval.0 <= alignment.pos)
        && (coord.start + coord.confidence_interval.0 >= alignment.pos)
    {
        matches = true;
    }
    if (coord.start - coord.confidence_interval.0 <= alignment.end)
        && (coord.start + coord.confidence_interval.0 >= alignment.end)
    {
        matches = true;
    }
    matches
}

fn phase_is_consistent(alignments: &[FwdStrandSplitReadSegment], allow_unphased: bool) -> bool {
    if allow_unphased {
        return true;
    }
    let mut phaseset_tags = HashSet::new();
    let mut haplotype_tags = HashSet::new();
    for align in alignments.iter() {
        if let Some(haplotype_tag) = align.haplotype_tag {
            haplotype_tags.insert(haplotype_tag);
        }
        if let Some(phaseset_tag) = align.phaseset_tag {
            phaseset_tags.insert(phaseset_tag);
        }
    }

    phaseset_tags.len() <= 1 && phaseset_tags.len() <= 1
}
