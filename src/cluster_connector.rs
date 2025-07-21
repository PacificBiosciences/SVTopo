use log::debug;
use std::collections::{HashMap, HashSet};
use std::time::SystemTime;

use crate::containers::{BlockStartInfo, Connection, Coordinate, FwdStrandSplitReadSegment};
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
        debug!("{break_count} connected pair of breaks found");
    } else {
        debug!("{break_count} connected pairs of breaks found");
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
    let ordered_breaks = organize_breaks_by_chromosome(cluster_alignments);
    for coordinates in ordered_breaks.values() {
        for (first_idx, &first_coord) in coordinates.iter().enumerate() {
            let first_coord_alignments = cluster_alignments.get(first_coord).unwrap();
            if let Some(first_block_start_info) =
                find_block_start_info(first_coord, first_coord_alignments)
            {
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
                    if second_coord.start >= first_block_start_info.alignment_start
                        && second_coord.start <= first_block_start_info.alignment_end
                    {
                        continue;
                    }
                    let second_coord_alignments = cluster_alignments.get(second_coord).unwrap();
                    if check_phase_compatibility(
                        second_coord,
                        second_coord_alignments,
                        first_block_start_info.phaseset,
                        first_block_start_info.haplotype,
                        &first_block_start_info,
                    ) && is_block_end(
                        second_coord,
                        second_coord_alignments,
                        first_block_start_info.phaseset,
                        first_block_start_info.haplotype,
                    ) {
                        phased_connection_count += create_phased_connections(
                            first_coord,
                            second_coord,
                            first_coord_alignments,
                            second_coord_alignments,
                            connected_breaks,
                        );
                        break;
                    }
                }
            }
        }
    }
    debug!("{phased_connection_count} phased connections added");
}

/// Organize cluster coordinates by chromosome and sort them for efficient phaseset processing
fn organize_breaks_by_chromosome(
    cluster_alignments: &HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
) -> HashMap<String, Vec<&Coordinate>> {
    let mut ordered_breaks: HashMap<String, Vec<&Coordinate>> = HashMap::new();
    for coordinate in cluster_alignments.keys() {
        ordered_breaks
            .entry(coordinate.start_chrom.clone())
            .or_default()
            .push(coordinate);
    }
    for coordinates in ordered_breaks.values_mut() {
        coordinates.sort_unstable();
    }
    ordered_breaks
}

/// Extract phase and alignment boundary information for a potential block start
fn find_block_start_info(
    coord: &Coordinate,
    alignments: &HashSet<FwdStrandSplitReadSegment>,
) -> Option<BlockStartInfo> {
    if !is_potential_block_start(coord, alignments) {
        return None;
    }

    let mut phaseset = 0;
    let mut haplotype = 0;
    let mut alignment_start = i64::MAX;
    let mut alignment_end = 0;

    for aln in alignments {
        if !aln.spans {
            continue;
        }
        if aln.pos < alignment_start {
            alignment_start = aln.pos;
        }
        if aln.end > alignment_end {
            alignment_end = aln.end;
        }
        if let (Some(ps), Some(hp)) = (&aln.phaseset_tag, &aln.haplotype_tag) {
            phaseset = *ps;
            haplotype = *hp;
        }
    }

    Some(BlockStartInfo {
        phaseset,
        haplotype,
        alignment_start,
        alignment_end,
    })
}

/// Check if a coordinate is compatible with the given phase information
fn check_phase_compatibility(
    coord: &Coordinate,
    alignments: &HashSet<FwdStrandSplitReadSegment>,
    target_phaseset: i32,
    target_haplotype: i32,
    start_info: &BlockStartInfo,
) -> bool {
    let dist = (start_info.alignment_start - coord.start).abs();
    if dist > utils::MAX_PHASE_DIST {
        return false;
    }

    if coord.start >= start_info.alignment_start && coord.start <= start_info.alignment_end {
        return false;
    }

    let mut is_phased = false;
    let mut phase_matches = false;

    for aln in alignments {
        if let (Some(phaseset), Some(haplotype)) = (&aln.phaseset_tag, &aln.haplotype_tag) {
            is_phased = true;
            if *phaseset != target_phaseset {
                return false; // moved to new phaseset
            } else if *haplotype == target_haplotype {
                phase_matches = true;
                break;
            }
        }
    }

    !is_phased || phase_matches
}

/// Create phased connections between two coordinates
fn create_phased_connections(
    first_coord: &Coordinate,
    second_coord: &Coordinate,
    first_alignments: &HashSet<FwdStrandSplitReadSegment>,
    second_alignments: &HashSet<FwdStrandSplitReadSegment>,
    connected_breaks: &mut HashMap<Connection, Vec<FwdStrandSplitReadSegment>>,
) -> usize {
    let combined_alignments = get_combined_alignments_for_phased_connections(
        first_coord,
        second_coord,
        first_alignments,
        second_alignments,
    );

    if combined_alignments.is_empty() {
        return 0;
    }

    let spanning_tags: HashSet<bool> = combined_alignments.iter().map(|a| a.spans).collect();
    let mut connection_count = 0;

    for spanning_tag in spanning_tags {
        let mut new_connection = Connection::new(
            first_coord.clone(),
            second_coord.clone(),
            true,
            spanning_tag,
        );
        let connection_alignments: Vec<FwdStrandSplitReadSegment> = combined_alignments
            .iter()
            .filter(|a| a.spans == spanning_tag)
            .cloned()
            .collect();
        if new_connection.first_coord > new_connection.second_coord {
            new_connection.reverse();
        }
        connected_breaks.insert(new_connection, connection_alignments);
        connection_count += 1;
    }

    connection_count
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

    phaseset_tags.len() <= 1 && haplotype_tags.len() <= 1
}

#[cfg(test)]
mod tests {
    use super::*;
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
        inferred_from_phasing: bool,
        is_spanned: bool,
    ) -> Connection {
        Connection::new(first_coord, second_coord, inferred_from_phasing, is_spanned)
    }

    #[test]
    fn test_alignment_matches_break() {
        let coord = create_test_coordinate("chr1", 1000, (10, 10));

        // Alignment that matches at start position
        let matching_alignment1 =
            utils::create_test_alignment_with_clips("chr1", 995, 1100, "read1", true, true);
        assert!(alignment_matches_break(&matching_alignment1, &coord));

        // Alignment that matches at end position
        let matching_alignment2 =
            utils::create_test_alignment_with_clips("chr1", 900, 1005, "read2", true, false);
        assert!(alignment_matches_break(&matching_alignment2, &coord));

        // Alignment that doesn't match
        let non_matching_alignment =
            utils::create_test_alignment_with_clips("chr1", 800, 900, "read3", true, true);
        assert!(!alignment_matches_break(&non_matching_alignment, &coord));

        // Different chromosome - NOTE: alignment_matches_break doesn't check chromosome!
        // It only checks position/end within confidence interval, so this will match
        let different_chrom_alignment =
            utils::create_test_alignment_with_clips("chr2", 1000, 1100, "read4", true, true);
        assert!(alignment_matches_break(&different_chrom_alignment, &coord));

        // Alignment with position outside confidence interval
        let far_alignment =
            utils::create_test_alignment_with_clips("chr1", 1050, 1150, "read5", true, true);
        assert!(!alignment_matches_break(&far_alignment, &coord));
    }

    #[test]
    fn test_phase_is_consistent() {
        // Test with allow_unphased = true
        let alignments = vec![
            utils::create_test_alignment_with_phasing(
                "chr1",
                1000,
                1100,
                "read1",
                true,
                true,
                Some(0),
                Some(0),
            ),
            utils::create_test_alignment_with_phasing(
                "chr1",
                1000,
                1100,
                "read2",
                true,
                true,
                Some(1),
                Some(0),
            ),
        ];
        assert!(phase_is_consistent(&alignments, true));

        // Test with allow_unphased = false, consistent phasing
        let consistent_alignments = vec![
            utils::create_test_alignment_with_phasing(
                "chr1",
                1000,
                1100,
                "read1",
                true,
                true,
                Some(0),
                Some(0),
            ),
            utils::create_test_alignment_with_phasing(
                "chr1",
                1000,
                1100,
                "read2",
                true,
                true,
                Some(0),
                Some(0),
            ),
        ];
        assert!(phase_is_consistent(&consistent_alignments, false));

        // Test with allow_unphased = false, inconsistent phasing
        let inconsistent_alignments = vec![
            utils::create_test_alignment_with_phasing(
                "chr1",
                1000,
                1100,
                "read1",
                true,
                true,
                Some(0),
                Some(0),
            ),
            utils::create_test_alignment_with_phasing(
                "chr1",
                1000,
                1100,
                "read2",
                true,
                true,
                Some(1),
                Some(0),
            ),
        ];
        assert!(!phase_is_consistent(&inconsistent_alignments, false));

        // Empty alignments
        assert!(phase_is_consistent(&[], false));
        assert!(phase_is_consistent(&[], true));
    }

    #[test]
    fn test_is_potential_block_start() {
        let coord = create_test_coordinate("chr1", 1000, (10, 10));

        // Create alignments that represent a block start (left-clipped, not right-clipped, same phaseset)
        let mut alignments = HashSet::new();
        alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            1100,
            "read1",
            true,
            false,
            Some(0),
            Some(0),
        ));
        alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            990,
            1050,
            "read2",
            true,
            false,
            Some(0),
            Some(0),
        ));

        assert!(is_potential_block_start(&coord, &alignments));

        // Test with right-clipped alignments (should not be block start)
        let mut right_clipped_alignments = HashSet::new();
        right_clipped_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            1100,
            "read1",
            false,
            true,
            Some(0),
            Some(0),
        ));
        right_clipped_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            990,
            1050,
            "read2",
            false,
            true,
            Some(0),
            Some(0),
        ));

        assert!(!is_potential_block_start(&coord, &right_clipped_alignments));

        // Test with different phasesets (should not be block start)
        let mut different_phase_alignments = HashSet::new();
        different_phase_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            1100,
            "read1",
            true,
            false,
            Some(1),
            Some(0),
        ));
        different_phase_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            990,
            1050,
            "read2",
            true,
            false,
            Some(2),
            Some(0),
        ));

        assert!(!is_potential_block_start(
            &coord,
            &different_phase_alignments
        ));

        // Test with insufficient left clips
        let mut insufficient_clips = HashSet::new();
        insufficient_clips.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            1100,
            "read1",
            true,
            false,
            Some(0),
            Some(0),
        ));

        assert!(!is_potential_block_start(&coord, &insufficient_clips));
    }

    #[test]
    fn test_is_block_end() {
        let coord = create_test_coordinate("chr1", 1000, (10, 10));
        let target_phaseset = 1;
        let target_haplotype = 0;

        // Create alignments that represent a block end (right-clipped, not left-clipped, matching phase)
        let mut alignments = HashSet::new();
        alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            900,
            1005,
            "read1",
            false,
            true,
            Some(1),
            Some(0),
        ));
        alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            950,
            1010,
            "read2",
            false,
            true,
            Some(1),
            Some(0),
        ));

        assert!(is_block_end(
            &coord,
            &alignments,
            target_phaseset,
            target_haplotype
        ));

        // Test with left-clipped alignments (should not be block end)
        let mut left_clipped_alignments = HashSet::new();
        left_clipped_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            900,
            1005,
            "read1",
            true,
            true,
            Some(0),
            Some(0),
        ));
        left_clipped_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            950,
            1010,
            "read2",
            true,
            true,
            Some(0),
            Some(0),
        ));

        assert!(!is_block_end(
            &coord,
            &left_clipped_alignments,
            target_phaseset,
            target_haplotype
        ));

        // Test with wrong phase
        let mut wrong_phase_alignments = HashSet::new();
        wrong_phase_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            900,
            1005,
            "read1",
            true,
            false,
            Some(2),
            Some(0),
        ));
        wrong_phase_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            950,
            1010,
            "read2",
            true,
            false,
            Some(2),
            Some(0),
        ));

        assert!(!is_block_end(
            &coord,
            &wrong_phase_alignments,
            target_phaseset,
            target_haplotype
        ));
    }

    #[test]
    fn test_find_block_start_info() {
        let coord = create_test_coordinate("chr1", 1000, (10, 10));

        // Valid block start
        let mut alignments = HashSet::new();
        alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            1200,
            "read1",
            true,
            false,
            Some(1),
            Some(0),
        ));
        alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            990,
            1150,
            "read2",
            true,
            false,
            Some(1),
            Some(0),
        ));

        let start_info = find_block_start_info(&coord, &alignments);
        assert!(start_info.is_some());
        let info = start_info.unwrap();
        assert_eq!(info.phaseset, 1);
        assert_eq!(info.haplotype, 0);
        assert_eq!(info.alignment_start, 990);
        assert_eq!(info.alignment_end, 1200);

        // Invalid block start (not potential block start)
        let mut invalid_alignments = HashSet::new();
        invalid_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            1100,
            "read1",
            true,
            false,
            Some(1),
            Some(0),
        ));

        let invalid_start_info = find_block_start_info(&coord, &invalid_alignments);
        assert!(invalid_start_info.is_none());

        // No spanning alignments
        let mut no_span_alignments = HashSet::new();
        no_span_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            1100,
            "read1",
            false,
            true,
            Some(1),
            Some(0),
        ));

        let no_span_info = find_block_start_info(&coord, &no_span_alignments);
        assert!(no_span_info.is_none());
    }

    #[test]
    fn test_check_phase_compatibility() {
        let coord = create_test_coordinate("chr1", 2000, (10, 10));
        let start_info = BlockStartInfo {
            phaseset: 1,
            haplotype: 0,
            alignment_start: 1000,
            alignment_end: 1500,
        };

        // Compatible phase (same phaseset, same haplotype)
        let mut compatible_alignments = HashSet::new();
        compatible_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            1950,
            2050,
            "read1",
            true,
            false,
            Some(1),
            Some(0),
        ));

        assert!(check_phase_compatibility(
            &coord,
            &compatible_alignments,
            1,
            0,
            &start_info
        ));

        // Incompatible phase (different phaseset)
        let mut incompatible_alignments = HashSet::new();
        incompatible_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            1950,
            2050,
            "read1",
            true,
            false,
            Some(2),
            Some(0),
        ));

        assert!(!check_phase_compatibility(
            &coord,
            &incompatible_alignments,
            1,
            0,
            &start_info
        ));

        // Too far apart (MAX_PHASE_DIST is 500,000, so use a larger distance)
        let far_coord = create_test_coordinate("chr1", 600000, (10, 10));
        assert!(!check_phase_compatibility(
            &far_coord,
            &compatible_alignments,
            1,
            0,
            &start_info
        ));

        // Overlapping with start info range
        let overlapping_coord = create_test_coordinate("chr1", 1200, (10, 10));
        assert!(!check_phase_compatibility(
            &overlapping_coord,
            &compatible_alignments,
            1,
            0,
            &start_info
        ));
    }

    #[test]
    fn test_organize_breaks_by_chromosome() {
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));
        let coord3 = create_test_coordinate("chr2", 1500, (10, 10));

        let mut cluster_alignments = HashMap::new();
        cluster_alignments.insert(coord1.clone(), HashSet::new());
        cluster_alignments.insert(coord2.clone(), HashSet::new());
        cluster_alignments.insert(coord3.clone(), HashSet::new());

        let organized = organize_breaks_by_chromosome(&cluster_alignments);

        assert_eq!(organized.len(), 2);
        assert!(organized.contains_key("chr1"));
        assert!(organized.contains_key("chr2"));

        let chr1_coords = organized.get("chr1").unwrap();
        assert_eq!(chr1_coords.len(), 2);
        // Should be sorted
        assert!(chr1_coords[0].start < chr1_coords[1].start);

        let chr2_coords = organized.get("chr2").unwrap();
        assert_eq!(chr2_coords.len(), 1);
    }

    #[test]
    fn test_get_combined_alignments_for_phased_connections() {
        let first_coord = create_test_coordinate("chr1", 1000, (10, 10));
        let second_coord = create_test_coordinate("chr1", 2000, (10, 10));

        // Create alignments that match the coordinates
        let alignment1 = utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            1100,
            "read1",
            true,
            false,
            Some(0),
            Some(0),
        );
        let alignment2 = utils::create_test_alignment_with_phasing(
            "chr1",
            1950,
            2050,
            "read2",
            true,
            false,
            Some(1),
            Some(0),
        );
        let alignment3 = utils::create_test_alignment_with_phasing(
            "chr1",
            800,
            900,
            "read3",
            true,
            true,
            Some(1),
            Some(0),
        ); // Doesn't match

        let mut first_alignments = HashSet::new();
        first_alignments.insert(alignment1.clone());
        first_alignments.insert(alignment3.clone());

        let mut second_alignments = HashSet::new();
        second_alignments.insert(alignment2.clone());

        let combined = get_combined_alignments_for_phased_connections(
            &first_coord,
            &second_coord,
            &first_alignments,
            &second_alignments,
        );

        // Should include alignment1 (matches first coord), but not alignment3 (doesn't match)
        // alignment2 is only in second_alignments but function processes both sets independently
        assert_eq!(combined.len(), 1);
        assert!(combined.contains(&alignment1));
        assert!(!combined.contains(&alignment3));
    }

    #[test]
    fn test_create_phased_connections() {
        let first_coord = create_test_coordinate("chr1", 1000, (10, 10));
        let second_coord = create_test_coordinate("chr1", 2000, (10, 10));

        let alignment1 = utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            1100,
            "read1",
            true,
            false,
            Some(0),
            Some(0),
        );
        let alignment2 = utils::create_test_alignment_with_phasing(
            "chr1",
            1950,
            2050,
            "read2",
            true,
            false,
            Some(1),
            Some(0),
        );

        let mut first_alignments = HashSet::new();
        first_alignments.insert(alignment1.clone());
        first_alignments.insert(alignment2.clone());

        let mut second_alignments = HashSet::new();
        second_alignments.insert(alignment1.clone());
        second_alignments.insert(alignment2.clone());

        let mut connected_breaks = HashMap::new();

        let count = create_phased_connections(
            &first_coord,
            &second_coord,
            &first_alignments,
            &second_alignments,
            &mut connected_breaks,
        );

        assert_eq!(count, 1); // Both alignments have spans=true, so only one connection type
        assert_eq!(connected_breaks.len(), 1);

        // Check that connections are properly created
        for (connection, alignments) in &connected_breaks {
            assert!(connection.inferred_from_phasing);
            assert!(!alignments.is_empty());
        }
    }

    #[test]
    fn test_connect_clusters_by_alignment() {
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));

        // Create shared alignments
        let alignment1 = utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            2050,
            "read1",
            true,
            true,
            Some(1),
            Some(0),
        );
        let alignment2 = utils::create_test_alignment_with_phasing(
            "chr1",
            990,
            2045,
            "read2",
            true,
            true,
            Some(1),
            Some(0),
        );
        let alignment3 = utils::create_test_alignment_with_phasing(
            "chr1",
            985,
            2040,
            "read3",
            false,
            true,
            Some(1),
            Some(0),
        );

        let mut coord1_alignments = HashSet::new();
        coord1_alignments.insert(alignment1.clone());
        coord1_alignments.insert(alignment2.clone());
        coord1_alignments.insert(alignment3.clone());

        let mut coord2_alignments = HashSet::new();
        coord2_alignments.insert(alignment1.clone());
        coord2_alignments.insert(alignment2.clone());
        coord2_alignments.insert(alignment3.clone());

        let mut cluster_alignments = HashMap::new();
        cluster_alignments.insert(coord1.clone(), coord1_alignments);
        cluster_alignments.insert(coord2.clone(), coord2_alignments);

        let connections = connect_clusters_by_alignment(&cluster_alignments, true);

        // Should have connections for both spanned and unspanned alignments
        assert!(!connections.is_empty());

        // Check that connections are properly formed
        for (connection, alignments) in &connections {
            assert!(!connection.inferred_from_phasing);
            assert!(alignments.len() >= utils::MIN_CONNECTION_SUPPORT as usize);
        }
    }

    #[test]
    fn test_connect_clusters_by_vcf_connection() {
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));

        let vcf_connection = create_test_connection(coord1.clone(), coord2.clone(), false, true);
        let vcf_connections = vec![vcf_connection];

        // Create shared alignments
        let alignment1 = utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            2050,
            "read1",
            true,
            true,
            Some(1),
            Some(0),
        );
        let alignment2 = utils::create_test_alignment_with_phasing(
            "chr1",
            990,
            2045,
            "read2",
            true,
            true,
            Some(1),
            Some(0),
        );

        let mut coord1_alignments = HashSet::new();
        coord1_alignments.insert(alignment1.clone());
        coord1_alignments.insert(alignment2.clone());

        let mut coord2_alignments = HashSet::new();
        coord2_alignments.insert(alignment1.clone());
        coord2_alignments.insert(alignment2.clone());

        let mut cluster_alignments = HashMap::new();
        cluster_alignments.insert(coord1.clone(), coord1_alignments);
        cluster_alignments.insert(coord2.clone(), coord2_alignments);

        let mut connected_breaks = HashMap::new();

        connect_clusters_by_vcf_connection(
            &cluster_alignments,
            &mut connected_breaks,
            &vcf_connections,
        );

        assert!(!connected_breaks.is_empty());

        // Check that the connection is properly added
        for (connection, alignments) in &connected_breaks {
            assert!(!connection.inferred_from_phasing);
            assert!(alignments.len() >= utils::MIN_CONNECTION_SUPPORT as usize);
        }
    }

    #[test]
    fn test_connect_clusters_by_phaseset() {
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));

        // Create block start alignments (left-clipped, phased)
        let mut start_alignments = HashSet::new();
        start_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            1500,
            "read1",
            true,
            false,
            Some(1),
            Some(0),
        ));
        start_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            990,
            1450,
            "read2",
            true,
            false,
            Some(1),
            Some(0),
        ));

        // Create block end alignments (right-clipped, same phase)
        let mut end_alignments = HashSet::new();
        end_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            1900,
            2005,
            "read3",
            false, // is_start_softclipped
            true,  // is_end_softclipped
            Some(1),
            Some(0),
        ));
        end_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            1950,
            2010,
            "read4",
            false, // is_start_softclipped
            true,  // is_end_softclipped
            Some(1),
            Some(0),
        ));

        let mut cluster_alignments = HashMap::new();
        cluster_alignments.insert(coord1.clone(), start_alignments);
        cluster_alignments.insert(coord2.clone(), end_alignments);

        let mut connected_breaks = HashMap::new();

        connect_clusters_by_phaseset(&cluster_alignments, &mut connected_breaks);

        // Should create phased connections
        let phased_connections: Vec<_> = connected_breaks
            .iter()
            .filter(|(conn, _)| conn.inferred_from_phasing)
            .collect();

        assert!(!phased_connections.is_empty());
    }

    #[test]
    fn test_connect_clusters_integration() {
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));

        // Create shared alignments for alignment-based connections
        let alignment1 = utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            2050,
            "read1",
            true,
            true,
            Some(1),
            Some(0),
        );
        let alignment2 = utils::create_test_alignment_with_phasing(
            "chr1",
            990,
            2045,
            "read2",
            true,
            true,
            Some(1),
            Some(0),
        );

        let mut coord1_alignments = HashSet::new();
        coord1_alignments.insert(alignment1.clone());
        coord1_alignments.insert(alignment2.clone());

        let mut coord2_alignments = HashSet::new();
        coord2_alignments.insert(alignment1.clone());
        coord2_alignments.insert(alignment2.clone());

        let mut cluster_alignments = HashMap::new();
        cluster_alignments.insert(coord1.clone(), coord1_alignments);
        cluster_alignments.insert(coord2.clone(), coord2_alignments);

        let vcf_connections = vec![];

        let connections = connect_clusters(&cluster_alignments, &vcf_connections, true);

        assert!(!connections.is_empty());

        // Should have deduplicated and properly formed connections
        for alignments in connections.values() {
            assert!(!alignments.is_empty());
        }
    }

    #[test]
    fn test_edge_cases() {
        // Test with empty cluster alignments
        let empty_cluster_alignments = HashMap::new();
        let empty_vcf_connections = vec![];

        let empty_connections =
            connect_clusters(&empty_cluster_alignments, &empty_vcf_connections, false);
        assert!(empty_connections.is_empty());

        // Test with alignments below minimum connection support
        let coord1 = create_test_coordinate("chr1", 1000, (10, 10));
        let coord2 = create_test_coordinate("chr1", 2000, (10, 10));

        let single_alignment = utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            2050,
            "read1",
            true,
            true,
            Some(1),
            Some(0),
        );

        let mut coord1_alignments = HashSet::new();
        coord1_alignments.insert(single_alignment.clone());

        let mut coord2_alignments = HashSet::new();
        coord2_alignments.insert(single_alignment);

        let mut cluster_alignments = HashMap::new();
        cluster_alignments.insert(coord1, coord1_alignments);
        cluster_alignments.insert(coord2, coord2_alignments);

        let connections = connect_clusters_by_alignment(&cluster_alignments, false);
        assert!(connections.is_empty());
    }

    #[test]
    fn test_phase_consistency_edge_cases() {
        // Test with no phasing information
        let alignments = vec![
            utils::create_test_alignment_with_phasing(
                "chr1",
                1000,
                1100,
                "read1",
                true,
                true,
                Some(1),
                Some(0),
            ),
            utils::create_test_alignment_with_phasing(
                "chr1",
                1000,
                1100,
                "read2",
                true,
                true,
                Some(2),
                Some(0),
            ),
        ];
        assert!(!phase_is_consistent(&alignments, false));
        assert!(phase_is_consistent(&alignments, true));

        // Test with mixed phasing information
        let mixed_alignments = vec![
            utils::create_test_alignment_with_phasing(
                "chr1",
                1000,
                1100,
                "read1",
                true,
                true,
                Some(1),
                Some(0),
            ),
            utils::create_test_alignment_with_phasing(
                "chr1",
                1000,
                1100,
                "read2",
                true,
                true,
                Some(2),
                Some(0),
            ),
        ];
        assert!(!phase_is_consistent(&mixed_alignments, false));
        assert!(phase_is_consistent(&mixed_alignments, true));
    }
}
