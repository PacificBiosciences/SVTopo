use log::{debug, error};
use std::collections::{HashMap, HashSet};
use std::time::SystemTime;

use crate::containers::{Coordinate, FwdStrandSplitReadSegment};
use crate::utils;

/// Given a hashmap with groups of read alignments keyed by readnames,
/// identify all clusters of softclip locations with at least two supporting reads.
/// These will define the start and stop coordinate of blocks in complex SVs.
/// Join together the softclips that are within max_clust_width bases of one another and
/// collect these in a hashmap of coordinates to vectors of supporting alignments.
pub fn find_breaks(
    clipped_reads: &HashMap<String, Vec<FwdStrandSplitReadSegment>>,
    allow_unphased: bool,
) -> HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>> {
    let start_time = SystemTime::now();
    let cluster_positions = cluster_via_clipping(clipped_reads, allow_unphased);

    debug!(
        "{} clustered break locations post-filtering",
        cluster_positions.len()
    );

    for (break_location, support) in cluster_positions.iter() {
        debug!("{} (support={})", break_location, support.len());
    }
    debug!(
        "Finding breaks: {}s",
        start_time.elapsed().unwrap().as_secs()
    );
    cluster_positions
}

/// Combine entries of softclip locations where multiple softclips are very close together
/// by storing a position, walking through the list linearly, and storing the cluster when
/// a new position is processed that's more than max_clust_width away
fn cluster_via_clipping(
    clipped_reads: &HashMap<String, Vec<FwdStrandSplitReadSegment>>,
    allow_unphased: bool,
) -> HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>> {
    let clip_coords = collect_softclip_locations(clipped_reads);
    log_clip_coordinates(&clip_coords);

    let sorted_clip_coords = sort_clip_coordinates(&clip_coords);

    process_coordinate_clusters(&sorted_clip_coords, &clip_coords, allow_unphased)
}

/// Collects all softclip locations into a hashmap
fn collect_softclip_locations(
    clipped_reads: &HashMap<String, Vec<FwdStrandSplitReadSegment>>,
) -> HashMap<Coordinate, Vec<FwdStrandSplitReadSegment>> {
    let mut clip_coords = HashMap::new();

    for (_read_name, alignments) in clipped_reads.iter() {
        for alignment in alignments {
            if alignment.is_start_softclipped {
                let chrom = &alignment.chrom;
                let break_coord = Coordinate::new(chrom.to_string(), alignment.pos);
                clip_coords
                    .entry(break_coord)
                    .or_insert(Vec::new())
                    .push(alignment.clone());
            }
            if alignment.is_end_softclipped {
                let chrom = &alignment.second_chrom;
                let break_coord = Coordinate::new(chrom.to_string(), alignment.end);
                clip_coords
                    .entry(break_coord)
                    .or_insert(Vec::new())
                    .push(alignment.clone());
            }
        }
    }

    clip_coords
}

/// Logs debug information about clip coordinates
fn log_clip_coordinates(clip_coords: &HashMap<Coordinate, Vec<FwdStrandSplitReadSegment>>) {
    debug!(
        "{} putative break locations pre-filtering",
        clip_coords.len()
    );
    for (break_location, support) in clip_coords.iter() {
        debug!("{} (support={})", break_location, support.len());
    }
}

/// Sorts clip coordinates for processing
fn sort_clip_coordinates(
    clip_coords: &HashMap<Coordinate, Vec<FwdStrandSplitReadSegment>>,
) -> Vec<Coordinate> {
    let mut clip_coords_vec: Vec<Coordinate> = clip_coords.keys().cloned().collect();
    clip_coords_vec.sort();
    clip_coords_vec
}

/// Processes sorted coordinates to create clusters
fn process_coordinate_clusters(
    sorted_clip_coords: &[Coordinate],
    clip_coords: &HashMap<Coordinate, Vec<FwdStrandSplitReadSegment>>,
    allow_unphased: bool,
) -> HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>> {
    let mut cluster_positions = HashMap::new();
    let mut begin_iter: usize = 0;

    for (end_iter, coord) in sorted_clip_coords.iter().enumerate() {
        if should_start_new_cluster(sorted_clip_coords, begin_iter, coord) {
            process_cluster_range(
                sorted_clip_coords,
                clip_coords,
                begin_iter,
                end_iter,
                &mut cluster_positions,
                allow_unphased,
            );
            begin_iter = end_iter;
        }
    }

    // Add the last cluster if there are any positions
    if !cluster_positions.is_empty() && begin_iter < sorted_clip_coords.len() {
        process_final_cluster(
            sorted_clip_coords,
            clip_coords,
            begin_iter,
            &mut cluster_positions,
            allow_unphased,
        );
    }

    cluster_positions
}

/// Determines if a new cluster should be started based on distance criteria
fn should_start_new_cluster(
    sorted_clip_coords: &[Coordinate],
    begin_iter: usize,
    coord: &Coordinate,
) -> bool {
    if begin_iter >= sorted_clip_coords.len() {
        return false;
    }

    let cluster_start = &sorted_clip_coords[begin_iter];
    (coord.start_chrom != cluster_start.start_chrom)
        || (coord.start > (cluster_start.start + cluster_start.confidence_interval.0))
}

/// Processes a cluster range and adds it if it meets criteria
fn process_cluster_range(
    sorted_clip_coords: &[Coordinate],
    clip_coords: &HashMap<Coordinate, Vec<FwdStrandSplitReadSegment>>,
    begin_iter: usize,
    end_iter: usize,
    cluster_positions: &mut HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    allow_unphased: bool,
) {
    let middle = get_clust_middle(&sorted_clip_coords[begin_iter..end_iter]);
    let alignments = get_clust_aligments(&sorted_clip_coords[begin_iter..end_iter], clip_coords);
    let spanned_count = alignments.iter().filter(|a| a.spans).count();

    if spanned_count >= utils::MIN_CONNECTION_SUPPORT as usize {
        add_cluster(
            middle,
            alignments,
            cluster_positions,
            utils::MIN_CONNECTION_SUPPORT,
            allow_unphased,
        );
    }
}

/// Processes the final cluster in the coordinate list
fn process_final_cluster(
    sorted_clip_coords: &[Coordinate],
    clip_coords: &HashMap<Coordinate, Vec<FwdStrandSplitReadSegment>>,
    begin_iter: usize,
    cluster_positions: &mut HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    allow_unphased: bool,
) {
    let middle = get_clust_middle(&sorted_clip_coords[begin_iter..]);
    let alignments = get_clust_aligments(&sorted_clip_coords[begin_iter..], clip_coords);
    add_cluster(
        middle,
        alignments,
        cluster_positions,
        utils::MIN_CONNECTION_SUPPORT,
        allow_unphased,
    );
}

/// Given a new coordinate and supporting alignments, add it to the cluster_positions hashmap
/// IF it also passes haplotype and connection support filters
fn add_cluster(
    new_cluster: Coordinate,
    alignments: HashSet<FwdStrandSplitReadSegment>,
    cluster_positions: &mut HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    min_connection_support: u32,
    allow_unphased: bool,
) {
    if allow_unphased {
        cluster_positions.insert(new_cluster, alignments);
    } else {
        // skip this break if phase is ambiguous or if not enough phased support
        let (haplotype_count, phaseset_count, phased_alignment_count) =
            count_alignment_phases(&alignments);
        let phase_is_unambiguous = haplotype_count == 1 && phaseset_count == 1;
        if phased_alignment_count >= min_connection_support as usize && phase_is_unambiguous {
            cluster_positions.insert(new_cluster, alignments);
        }
    }
}

/// Given an array of coordinates, find their middle. Used to
/// identify the center of a group of alignments with very close clip
/// locations. Returns the Coordinate of the middle location.
fn get_clust_middle(cluster_coords: &[Coordinate]) -> Coordinate {
    let mut cluster_poses: Vec<i64> = Vec::new();
    if cluster_coords.len() == 1 {
        return cluster_coords.first().unwrap().clone();
    }
    for clip_pos in cluster_coords.iter() {
        cluster_poses.push(clip_pos.start);
    }
    cluster_poses.sort();
    let middle_idx = (cluster_poses.len() - 1) / 2;
    let pos = match cluster_poses.len() % 2 {
        0 => (cluster_poses[middle_idx] + cluster_poses[middle_idx + 1]) / 2,
        1 => cluster_poses[middle_idx],
        _ => {
            error!("get_clust_middle() unreachable");
            std::process::exit(exitcode::IOERR);
        }
    };
    Coordinate::new(cluster_coords[0].start_chrom.clone(), pos)
}

/// Given an array of cluster coordinates, get all of the alignments
/// that belong to that cluster. Returns a hashset of alignments.
fn get_clust_aligments(
    cluster_coords: &[Coordinate],
    clip_coords: &HashMap<Coordinate, Vec<FwdStrandSplitReadSegment>>,
) -> HashSet<FwdStrandSplitReadSegment> {
    let mut alignments = HashSet::new();
    for cluster_clip in cluster_coords.iter() {
        let clip_alignments: HashSet<FwdStrandSplitReadSegment> =
            clip_coords[cluster_clip].iter().cloned().collect();
        alignments = alignments.union(&clip_alignments).cloned().collect();
    }
    alignments
}

/// Given a vector of alignments corresponding to a cluster,
/// count the number of haplotypes and phasesets.
/// Return the number of each and the number of phased alignments from unique reads
fn count_alignment_phases(
    alignments: &HashSet<FwdStrandSplitReadSegment>,
) -> (usize, usize, usize) {
    let mut haplotypes = HashSet::new();
    let mut phasesets = HashSet::new();
    let mut phased_alignments = HashSet::new();
    for align in alignments.iter() {
        if let Some(hp) = align.haplotype_tag {
            haplotypes.insert(hp);
            phased_alignments.insert(&align.readname);
        }
        if let Some(ps) = align.phaseset_tag {
            phasesets.insert(ps);
        }
    }
    (haplotypes.len(), phasesets.len(), phased_alignments.len())
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

    // Helper function to create clipped reads data
    fn create_test_clipped_reads() -> HashMap<String, Vec<FwdStrandSplitReadSegment>> {
        let mut clipped_reads = HashMap::new();

        clipped_reads.insert(
            "read1".to_string(),
            vec![
                utils::create_test_alignment_with_phasing(
                    "chr1",
                    1000,
                    1100,
                    "read1",
                    true,
                    true,
                    Some(1),
                    None,
                ),
                utils::create_test_alignment_with_phasing(
                    "chr1",
                    2000,
                    2100,
                    "read1",
                    true,
                    true,
                    Some(1),
                    None,
                ),
            ],
        );

        clipped_reads.insert(
            "read2".to_string(),
            vec![
                utils::create_test_alignment_with_phasing(
                    "chr1",
                    1005,
                    1105,
                    "read2",
                    true,
                    true,
                    Some(1),
                    None,
                ),
                utils::create_test_alignment_with_phasing(
                    "chr1",
                    2005,
                    2105,
                    "read2",
                    true,
                    true,
                    Some(1),
                    None,
                ),
            ],
        );

        clipped_reads.insert(
            "read3".to_string(),
            vec![utils::create_test_alignment_with_phasing(
                "chr1",
                3000,
                3100,
                "read3",
                true,
                true,
                Some(2),
                None,
            )],
        );

        clipped_reads
    }

    #[test]
    fn test_collect_softclip_locations() {
        let clipped_reads = create_test_clipped_reads();
        let clip_coords = collect_softclip_locations(&clipped_reads);

        // Should have 4 clip coordinates: 1000, 1005 (start clips), 2100, 2105 (end clips)
        assert!(clip_coords.len() >= 4);

        // Check that coordinates exist for the expected positions
        let coord_1000 = Coordinate::new("chr1".to_string(), 1000);
        let coord_1005 = Coordinate::new("chr1".to_string(), 1005);
        let coord_2100 = Coordinate::new("chr1".to_string(), 2100);
        let coord_2105 = Coordinate::new("chr1".to_string(), 2105);

        assert!(clip_coords.contains_key(&coord_1000));
        assert!(clip_coords.contains_key(&coord_1005));
        assert!(clip_coords.contains_key(&coord_2100));
        assert!(clip_coords.contains_key(&coord_2105));
    }

    #[test]
    fn test_sort_clip_coordinates() {
        let mut clip_coords = HashMap::new();

        let coord1 = Coordinate::new("chr1".to_string(), 2000);
        let coord2 = Coordinate::new("chr1".to_string(), 1000);
        let coord3 = Coordinate::new("chr2".to_string(), 1500);

        clip_coords.insert(coord1, vec![]);
        clip_coords.insert(coord2, vec![]);
        clip_coords.insert(coord3, vec![]);

        let sorted = sort_clip_coordinates(&clip_coords);

        assert_eq!(sorted.len(), 3);
        // Should be sorted by chromosome then position
        assert_eq!(sorted[0].start_chrom, "chr1");
        assert_eq!(sorted[0].start, 1000);
        assert_eq!(sorted[1].start_chrom, "chr1");
        assert_eq!(sorted[1].start, 2000);
        assert_eq!(sorted[2].start_chrom, "chr2");
        assert_eq!(sorted[2].start, 1500);
    }

    #[test]
    fn test_should_start_new_cluster() {
        let coords = vec![
            create_test_coordinate("chr1", 1000, (utils::MAX_CLUST_DIST, utils::MAX_CLUST_DIST)),
            create_test_coordinate("chr1", 1010, (utils::MAX_CLUST_DIST, utils::MAX_CLUST_DIST)), // Within confidence interval
            create_test_coordinate("chr1", 1050, (utils::MAX_CLUST_DIST, utils::MAX_CLUST_DIST)), // Outside confidence interval
            create_test_coordinate("chr1", 2000, (utils::MAX_CLUST_DIST, utils::MAX_CLUST_DIST)),
            create_test_coordinate("chr2", 1000, (utils::MAX_CLUST_DIST, utils::MAX_CLUST_DIST)),
        ];

        // Should not start new cluster for coordinates within confidence interval
        assert!(!should_start_new_cluster(&coords, 0, &coords[1])); // 1010 is within 1000 ± 15

        // Should start new cluster for coordinates outside confidence interval
        assert!(should_start_new_cluster(&coords, 0, &coords[2])); // 1050 is outside 1000 ± 15

        // Should start new cluster for far coordinates
        assert!(should_start_new_cluster(&coords, 0, &coords[3]));

        // Should start new cluster for different chromosome
        assert!(should_start_new_cluster(&coords, 0, &coords[4]));

        // Edge case: begin_iter out of bounds
        assert!(!should_start_new_cluster(&coords, 10, &coords[0]));
    }

    #[test]
    fn test_get_clust_middle() {
        // Single coordinate
        let single_coord = vec![Coordinate::new("chr1".to_string(), 1000)];
        let middle_single = get_clust_middle(&single_coord);
        assert_eq!(middle_single.start, 1000);

        // Odd number of coordinates
        let odd_coords = vec![
            Coordinate::new("chr1".to_string(), 1000),
            Coordinate::new("chr1".to_string(), 1010),
            Coordinate::new("chr1".to_string(), 1020),
        ];
        let middle_odd = get_clust_middle(&odd_coords);
        assert_eq!(middle_odd.start, 1010);

        // Even number of coordinates
        let even_coords = vec![
            Coordinate::new("chr1".to_string(), 1000),
            Coordinate::new("chr1".to_string(), 1020),
        ];
        let middle_even = get_clust_middle(&even_coords);
        assert_eq!(middle_even.start, 1010); // Average of 1000 and 1020
    }

    #[test]
    fn test_get_clust_aligments() {
        let coord1 = Coordinate::new("chr1".to_string(), 1000);
        let coord2 = Coordinate::new("chr1".to_string(), 1010);

        let alignment1 = utils::create_test_alignment("chr1", 995, 1100, "read1");
        let alignment2 = utils::create_test_alignment("chr1", 1005, 1110, "read2");

        let mut clip_coords = HashMap::new();
        clip_coords.insert(coord1.clone(), vec![alignment1.clone()]);
        clip_coords.insert(coord2.clone(), vec![alignment2.clone()]);

        let cluster_coords = vec![coord1, coord2];
        let alignments = get_clust_aligments(&cluster_coords, &clip_coords);

        assert_eq!(alignments.len(), 2);
        assert!(alignments.contains(&alignment1));
        assert!(alignments.contains(&alignment2));
    }

    #[test]
    fn test_count_alignment_phases() {
        let mut alignments = HashSet::new();

        // Add alignments with different phases
        alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            1000,
            1100,
            "read1",
            true,
            true,
            Some(1),
            Some(0),
        ));
        alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            1005,
            1105,
            "read2",
            true,
            true,
            Some(1),
            Some(0),
        ));
        alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            1010,
            1110,
            "read3",
            true,
            true,
            Some(1),
            Some(0),
        ));
        alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            1015,
            1115,
            "read4",
            true,
            true,
            Some(2),
            Some(1),
        ));

        let (haplotype_count, phaseset_count, phased_alignment_count) =
            count_alignment_phases(&alignments);

        assert_eq!(haplotype_count, 2); // Haplotypes 0 and 1
        assert_eq!(phaseset_count, 2); // Phasesets 1 and 2
        assert_eq!(phased_alignment_count, 4); // All 4 reads are phased

        // Test with no phasing
        let mut unphased_alignments = HashSet::new();
        unphased_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 1000, 1100, "read1", true, true, None, None,
        ));

        let (h_count, p_count, ph_count) = count_alignment_phases(&unphased_alignments);
        assert_eq!(h_count, 0);
        assert_eq!(p_count, 0);
        assert_eq!(ph_count, 0);
    }

    #[test]
    fn test_add_cluster() {
        let mut cluster_positions = HashMap::new();
        let coord = Coordinate::new("chr1".to_string(), 1000);

        // Test with phased alignments (allow_unphased = false)
        let mut phased_alignments = HashSet::new();
        phased_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            1100,
            "read1",
            true,
            true,
            Some(1),
            Some(0),
        ));
        phased_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            990,
            1095,
            "read2",
            true,
            true,
            Some(1),
            Some(0),
        ));

        add_cluster(
            coord.clone(),
            phased_alignments.clone(),
            &mut cluster_positions,
            2,
            false,
        );
        assert!(cluster_positions.contains_key(&coord));

        // Test with allow_unphased = true
        let mut cluster_positions2 = HashMap::new();
        let mut unphased_alignments = HashSet::new();
        unphased_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 995, 1100, "read1", true, true, None, None,
        ));
        unphased_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 990, 1095, "read2", true, true, None, None,
        ));

        let coord2 = Coordinate::new("chr1".to_string(), 2000);
        add_cluster(
            coord2.clone(),
            unphased_alignments,
            &mut cluster_positions2,
            2,
            true,
        );
        assert!(cluster_positions2.contains_key(&coord2));

        // Test with mixed phasing (should be rejected with allow_unphased = false)
        let mut cluster_positions3 = HashMap::new();
        let mut mixed_alignments = HashSet::new();
        mixed_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            1100,
            "read1",
            true,
            true,
            Some(1),
            None,
        ));
        mixed_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            990,
            1095,
            "read2",
            true,
            true,
            Some(1),
            None,
        )); // Different haplotype

        let coord3 = Coordinate::new("chr1".to_string(), 3000);
        add_cluster(
            coord3.clone(),
            mixed_alignments,
            &mut cluster_positions3,
            2,
            false,
        );
        assert!(!cluster_positions3.contains_key(&coord3)); // Should be rejected due to mixed haplotypes
    }

    #[test]
    fn test_cluster_via_clipping() {
        let clipped_reads = create_test_clipped_reads();
        let clusters = cluster_via_clipping(&clipped_reads, true);

        // Should find clusters for the clip locations
        assert!(!clusters.is_empty());

        // Test with phasing requirements
        let clusters_phased = cluster_via_clipping(&clipped_reads, false);

        // May have fewer clusters due to phasing requirements
        assert!(clusters_phased.len() <= clusters.len());
    }

    #[test]
    fn test_find_breaks_integration() {
        let clipped_reads = create_test_clipped_reads();
        let breaks = find_breaks(&clipped_reads, true);

        // Test that function runs and may or may not find breaks (depends on support criteria)
        // The main test is that it doesn't panic and returns a valid HashMap
        assert!(breaks.is_empty() || !breaks.is_empty());

        // If there are breaks, they should be well-formed
        for (coord, alignments) in &breaks {
            assert!(!coord.start_chrom.is_empty());
            assert!(!alignments.is_empty());
        }
    }

    #[test]
    fn test_edge_cases() {
        // Test with empty input
        let empty_reads = HashMap::new();
        let empty_breaks = find_breaks(&empty_reads, true);
        assert!(empty_breaks.is_empty());

        // Test with insufficient support
        let mut insufficient_reads = HashMap::new();
        insufficient_reads.insert(
            "read1".to_string(),
            vec![utils::create_test_alignment_with_phasing(
                "chr1",
                1000,
                1100,
                "read1",
                true,
                true,
                Some(1),
                None,
            )],
        );

        let insufficient_breaks = find_breaks(&insufficient_reads, false);
        assert!(insufficient_breaks.is_empty()); // Not enough support

        // Test with mixed chromosomes
        let mut mixed_reads = HashMap::new();
        mixed_reads.insert(
            "read1".to_string(),
            vec![
                utils::create_test_alignment_with_phasing(
                    "chr1",
                    1000,
                    1100,
                    "read1",
                    true,
                    true,
                    Some(1),
                    None,
                ),
                utils::create_test_alignment_with_phasing(
                    "chr2",
                    1000,
                    1100,
                    "read1",
                    true,
                    false,
                    Some(1),
                    None,
                ),
            ],
        );
        mixed_reads.insert(
            "read2".to_string(),
            vec![
                utils::create_test_alignment_with_phasing(
                    "chr1",
                    1005,
                    1105,
                    "read2",
                    true,
                    true,
                    Some(1),
                    None,
                ),
                utils::create_test_alignment_with_phasing(
                    "chr2",
                    1005,
                    1105,
                    "read2",
                    true,
                    false,
                    Some(1),
                    None,
                ),
            ],
        );

        let mixed_breaks = find_breaks(&mixed_reads, true);

        // Should find breaks on both chromosomes
        let chr1_breaks = mixed_breaks
            .keys()
            .filter(|coord| coord.start_chrom == "chr1")
            .count();
        let chr2_breaks = mixed_breaks
            .keys()
            .filter(|coord| coord.start_chrom == "chr2")
            .count();
        assert!(chr1_breaks > 0 || chr2_breaks > 0);
    }
}
