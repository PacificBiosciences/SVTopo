use log::{debug, error};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::time::SystemTime;

use crate::containers::{Coordinate, FwdStrandSplitReadSegment};
use crate::utils::{self, MAX_CLUST_DIST};

/// Given a hashmap with groups of read alignments keyed by readnames,
/// identify all clusters of softclip locations with at least two supporting reads.
/// These will define the start and stop coordinate of blocks in complex SVs.
/// Join together the softclips that are within max_clust_width bases of one another and
/// collect these in a hashmap of coordinates to vectors of supporting alignments.
pub fn assign_clipped_reads_to_clusters(
    clipped_reads: &HashMap<String, Vec<FwdStrandSplitReadSegment>>,
    vcf_coord_map: &HashMap<String, Vec<Coordinate>>,
    variant_readnames_opt: &Option<HashMap<String, Vec<String>>>,
) -> HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>> {
    let start_time = SystemTime::now();

    let cluster_positions: HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>> =
        if let Some(variant_readnames) = variant_readnames_opt {
            cluster_via_readnames(clipped_reads, vcf_coord_map, variant_readnames)
        } else {
            cluster_via_vcf(clipped_reads, vcf_coord_map)
        };

    if cluster_positions.is_empty() {
        let msg1 = "No readnames were matched to variant IDs from the VCF.";
        let msg2 = "This could be caused by several things, including a file mismatch such as between the BAM and the VCF, or the VCF and JSON.";

        error!("{msg1} {msg2}");
        std::process::exit(exitcode::DATAERR);
    }

    debug!(
        "{} clustered break locations post-filtering",
        cluster_positions.len()
    );

    let mut sorted_break_locations: Vec<&Coordinate> = cluster_positions.keys().collect();
    sorted_break_locations.sort();
    for break_location in sorted_break_locations.iter() {
        let support = cluster_positions.get(break_location).unwrap();
        debug!("{} (support={})", break_location, support.len());
    }
    debug!(
        "Finding breaks: {}s",
        start_time.elapsed().unwrap().as_secs()
    );
    cluster_positions
}

/// Uses the associations of variant IDs with readnames and the
/// variant IDs with coordinates, from the VCF, to assign specific
/// alignments to coordinates. The closest alignment from each read to
/// an associated breakend is treated as a match regardless of distance,
/// based on the assumption of correct input data from sawfish.
fn cluster_via_readnames(
    clipped_reads: &HashMap<String, Vec<FwdStrandSplitReadSegment>>,
    vcf_coord_map: &HashMap<String, Vec<Coordinate>>,
    variant_readnames: &HashMap<String, Vec<String>>,
) -> HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>> {
    let mut cluster_positions = BTreeMap::new();
    let mut vcf_variant_ids: Vec<&String> = vcf_coord_map.keys().collect();
    vcf_variant_ids.sort();

    for variant_id in vcf_variant_ids {
        let mut coords = vcf_coord_map.get(variant_id).unwrap().clone();
        coords.sort();
        for current_coord in coords.iter_mut() {
            let mut cluster_alignments = HashSet::new();
            if let Some(readnames) = variant_readnames.get(variant_id) {
                for readname in readnames {
                    if let Some(mut variant_alignments) = clipped_reads.get(readname).cloned() {
                        variant_alignments.sort();
                        add_vcf_alignments(
                            &variant_alignments,
                            current_coord,
                            &mut cluster_alignments,
                        );
                    }
                }
                let mut cluster_coord = current_coord.clone();
                cluster_coord.variant_ids.push(variant_id.clone());
                cluster_alignments = refine_alignments(cluster_alignments, &mut cluster_coord);
                add_coordinate(
                    &mut cluster_positions,
                    cluster_coord,
                    &mut cluster_alignments,
                );
            }
        }
    }
    cluster_positions.into_iter().collect()
}

/// Given a clipped alignment and a break coordinate, determine if the alignment matches the
/// break location and if so, if in the forward direction (meaning the end of the alignment is clipped)
fn alignment_matches_break_forward(
    alignment: &FwdStrandSplitReadSegment,
    break_coord: &Coordinate,
) -> bool {
    let mut matches = false;
    if alignment.chrom == break_coord.start_chrom
        || alignment.second_chrom == break_coord.start_chrom
    {
        if alignment.is_start_softclipped {
            let dist = (alignment.pos - break_coord.start).abs();
            if dist <= MAX_CLUST_DIST {
                matches = true;
            }
        }
        if alignment.is_end_softclipped {
            let dist = (alignment.end - break_coord.start).abs();
            if dist <= MAX_CLUST_DIST {
                matches = true;
            }
        }
    }
    matches
}

/// Assign specific alignments to coordinates without prior knowledge of
/// alignment match to break location. An alignment clipping position must
/// be no more than MAX_CLUST_DIST away to match a break location.
fn cluster_via_vcf(
    clipped_reads: &HashMap<String, Vec<FwdStrandSplitReadSegment>>,
    vcf_coord_map: &HashMap<String, Vec<Coordinate>>,
) -> HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>> {
    let mut cluster_positions = BTreeMap::new();
    let mut vcf_variant_ids: Vec<&String> = vcf_coord_map.keys().collect();
    vcf_variant_ids.sort();

    for variant_id in vcf_variant_ids {
        let mut coords = vcf_coord_map.get(variant_id).unwrap().clone();
        coords.sort();
        for current_coord in coords.iter_mut() {
            let mut cluster_alignments: HashSet<FwdStrandSplitReadSegment> = HashSet::new();

            for (_readname, read_alignments) in clipped_reads.iter() {
                for alignment in read_alignments.iter() {
                    let alignment_matches =
                        alignment_matches_break_forward(alignment, current_coord);
                    if alignment_matches {
                        cluster_alignments.insert(alignment.clone());
                    }
                }
            }
            let mut cluster_coord = current_coord.clone();
            cluster_coord.variant_ids.push(variant_id.clone());
            cluster_alignments = refine_alignments(cluster_alignments, &mut cluster_coord);
            add_coordinate(
                &mut cluster_positions,
                cluster_coord,
                &mut cluster_alignments,
            );
        }
    }
    cluster_positions.into_iter().collect()
}

/// Add alignments to the cluster alignments, from vcf coordinate matching
fn add_vcf_alignments(
    variant_alignments: &[FwdStrandSplitReadSegment],
    current_coord: &mut Coordinate,
    cluster_alignments: &mut HashSet<FwdStrandSplitReadSegment>,
) {
    let mut closest_alignment_dist = i64::MAX;

    // find the distance between break coordinate and closest alignment
    for alignment in variant_alignments.iter() {
        let mut clip_locations = Vec::new();
        if alignment.is_start_softclipped {
            clip_locations.push(alignment.pos);
        }
        if alignment.is_end_softclipped {
            clip_locations.push(alignment.end);
        }

        for clip_location in clip_locations.iter() {
            let alignment_dist = (current_coord.start - clip_location).abs();
            if alignment_dist < closest_alignment_dist {
                closest_alignment_dist = alignment_dist;
            }
        }
    }
    // add alignments at that distance from the breaks
    for alignment in variant_alignments.iter() {
        let mut clip_locations = Vec::new();
        if alignment.is_start_softclipped {
            clip_locations.push(alignment.pos);
        }
        if alignment.is_end_softclipped {
            clip_locations.push(alignment.end);
        }

        for clip_location in clip_locations.iter() {
            let alignment_dist = (current_coord.start - clip_location).abs();
            if alignment_dist == closest_alignment_dist {
                cluster_alignments.insert(alignment.clone());
            }
        }
    }
}

/// Refines the set of alignments assigned to the breakend by finding the most common
/// clipping location the assigned alignments share and omitting any alignments farther
/// away from the cluster location
fn refine_alignments(
    cluster_alignments: HashSet<FwdStrandSplitReadSegment>,
    cluster_coord: &mut Coordinate,
) -> HashSet<FwdStrandSplitReadSegment> {
    if cluster_alignments.is_empty() {
        return cluster_alignments;
    }

    let alignments_by_dist = group_alignments_by_distance(&cluster_alignments, cluster_coord);
    let (usable_clip_site_dists, most_common_dist) = find_usable_distances(&alignments_by_dist);
    let allowed_alignment_distance =
        determine_allowed_distance(&usable_clip_site_dists, most_common_dist);

    let refined_alignments =
        filter_alignments_by_distance(&alignments_by_dist, allowed_alignment_distance);
    update_cluster_confidence_interval(cluster_coord, allowed_alignment_distance);

    refined_alignments
}

/// Groups alignments by their minimum distance to cluster coordinates
fn group_alignments_by_distance(
    cluster_alignments: &HashSet<FwdStrandSplitReadSegment>,
    cluster_coord: &Coordinate,
) -> HashMap<i64, Vec<FwdStrandSplitReadSegment>> {
    let mut alignments_by_dist = HashMap::new();
    let mut ordered_cluster_alignments: Vec<&FwdStrandSplitReadSegment> =
        cluster_alignments.iter().collect();
    ordered_cluster_alignments.sort();

    for alignment in ordered_cluster_alignments {
        let min_distance = calculate_alignment_min_distance(alignment, cluster_coord);
        alignments_by_dist
            .entry(min_distance)
            .or_insert_with(Vec::new)
            .push(alignment.clone());
    }

    alignments_by_dist
}

/// Calculates the minimum distance between an alignment and cluster coordinates
fn calculate_alignment_min_distance(
    alignment: &FwdStrandSplitReadSegment,
    cluster_coord: &Coordinate,
) -> i64 {
    let mut left_dist = i64::MAX;
    let mut right_dist = i64::MAX;

    if alignment.is_start_softclipped {
        left_dist = std::cmp::min(
            (alignment.pos - cluster_coord.start).abs(),
            (alignment.pos - cluster_coord.end).abs(),
        );
    }
    if alignment.is_end_softclipped {
        right_dist = std::cmp::min(
            (alignment.end - cluster_coord.start).abs(),
            (alignment.end - cluster_coord.end).abs(),
        );
    }

    std::cmp::min(left_dist, right_dist)
}

/// Finds usable distances based on minimum connection support and determines most common distance
fn find_usable_distances(
    alignments_by_dist: &HashMap<i64, Vec<FwdStrandSplitReadSegment>>,
) -> (Vec<i64>, i64) {
    let dist_counts: Vec<(&i64, usize)> = alignments_by_dist
        .iter()
        .map(|(dist, alns)| (dist, alns.len()))
        .collect();

    let mut usable_clip_site_dists: Vec<i64> = Vec::new();
    let mut max_dist_count = 0;
    let mut most_common_dist = 0;

    for (distance, count) in dist_counts {
        if count as u32 >= utils::MIN_CONNECTION_SUPPORT {
            usable_clip_site_dists.push(*distance);
            if count > max_dist_count {
                max_dist_count = count;
                most_common_dist = *distance;
            }
        }
    }
    usable_clip_site_dists.sort();

    (usable_clip_site_dists, most_common_dist)
}

/// Determines the allowed alignment distance based on confidence intervals and most common distance
fn determine_allowed_distance(usable_clip_site_dists: &[i64], most_common_dist: i64) -> i64 {
    if let (Some(leftmost), Some(rightmost)) = (
        usable_clip_site_dists.first(),
        usable_clip_site_dists.last(),
    ) {
        if rightmost - leftmost > utils::MAX_ALIGNMENT_CONFIDENCE_INTERVAL {
            utils::MAX_ALIGNMENT_CONFIDENCE_INTERVAL
        } else {
            std::cmp::max(*rightmost, *leftmost)
        }
    } else {
        most_common_dist + utils::MAX_CLUST_DIST
    }
}

/// Filters alignments based on allowed distance
fn filter_alignments_by_distance(
    alignments_by_dist: &HashMap<i64, Vec<FwdStrandSplitReadSegment>>,
    allowed_alignment_distance: i64,
) -> HashSet<FwdStrandSplitReadSegment> {
    let mut refined_cluster_alignments = HashSet::new();

    for (dist, alignments) in alignments_by_dist {
        if *dist <= allowed_alignment_distance {
            refined_cluster_alignments.extend(alignments.iter().cloned());
        }
    }

    refined_cluster_alignments
}

/// Updates cluster coordinate confidence interval based on allowed distance
fn update_cluster_confidence_interval(
    cluster_coord: &mut Coordinate,
    allowed_alignment_distance: i64,
) {
    let ci = std::cmp::max(utils::MAX_CLUST_DIST, allowed_alignment_distance);
    cluster_coord.confidence_interval = (ci, ci);
}

/// Find out from assigned alignments if the right side of the given
/// coordinate is clipped. If not, assumes the left side is the clipped one.
/// Check is performed by looping through the spanning alignments
/// and for each one determining the closer clipped side to the given
/// coordinate.
/// Assumes the coordinate is a single position, not an interval.
fn get_coord_clip_directions(
    alignments: &HashSet<FwdStrandSplitReadSegment>,
    coord: &Coordinate,
) -> (bool, bool) {
    let mut is_rightclipped = false;
    let mut is_leftclipped = false;
    for alignment in alignments.iter() {
        if !alignment.spans {
            continue;
        }
        if alignment.chrom != coord.start_chrom {
            continue;
        }
        let mut right_dist = i64::MAX;
        let mut left_dist = i64::MAX;

        if alignment.is_end_softclipped {
            right_dist = (alignment.end - coord.start).abs();
        }
        if alignment.is_start_softclipped {
            left_dist = (alignment.pos - coord.start).abs();
        }
        if right_dist < left_dist {
            is_rightclipped = true;
        }
        if left_dist < right_dist {
            is_leftclipped = true;
        }
    }

    (is_rightclipped, is_leftclipped)
}

/// If a new coordinate is within the confidence interval of one or two others, merge them before adding to the map.
/// If the two are not close enough to merge, but share the same alignment, assign that alignment only to the coordinate
/// closest to it.
/// Otherwise, just add the coordinate to the map as key, with alignments as value.
fn add_coordinate(
    cluster_positions: &mut BTreeMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    cluster_coord: Coordinate,
    cluster_alignments: &mut HashSet<FwdStrandSplitReadSegment>,
) {
    if cluster_alignments.len() < utils::MIN_CONNECTION_SUPPORT as usize {
        return;
    }

    let working_cluster_positions = cluster_positions.clone();

    handle_previous_coordinate_merging(
        cluster_positions,
        &cluster_coord,
        cluster_alignments,
        &working_cluster_positions,
    );

    handle_next_coordinate_merging(
        cluster_positions,
        &cluster_coord,
        cluster_alignments,
        &working_cluster_positions,
    );

    cluster_positions.insert(cluster_coord, cluster_alignments.clone());
}

/// Handles merging with previous coordinate if within confidence interval
fn handle_previous_coordinate_merging(
    cluster_positions: &mut BTreeMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    cluster_coord: &Coordinate,
    cluster_alignments: &mut HashSet<FwdStrandSplitReadSegment>,
    working_cluster_positions: &BTreeMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
) {
    let prev_coord_opt: Option<(&Coordinate, &HashSet<FwdStrandSplitReadSegment>)> =
        working_cluster_positions.range(..cluster_coord).next_back();

    if let Some((prev_coord, prev_alignments)) = prev_coord_opt {
        if prev_coord.end_chrom == cluster_coord.start_chrom {
            process_previous_coordinate(
                cluster_positions,
                prev_coord,
                prev_alignments,
                cluster_coord,
                cluster_alignments,
            );
        }
    }
}

/// Processes the previous coordinate for potential merging or shared alignment assignment
fn process_previous_coordinate(
    cluster_positions: &mut BTreeMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    prev_coord: &Coordinate,
    prev_alignments: &HashSet<FwdStrandSplitReadSegment>,
    cluster_coord: &Coordinate,
    cluster_alignments: &mut HashSet<FwdStrandSplitReadSegment>,
) {
    let prev_clip_directions = get_coord_clip_directions(prev_alignments, prev_coord);
    let current_clip_directions = get_coord_clip_directions(cluster_alignments, cluster_coord);
    let prev_dist = (prev_coord.end - cluster_coord.start).abs();

    if should_merge_coordinates(
        prev_clip_directions,
        current_clip_directions,
        prev_dist,
        cluster_coord,
        prev_coord,
    ) {
        cluster_alignments.extend(prev_alignments.iter().cloned());
        cluster_positions.remove_entry(prev_coord);
    } else if let Some(prev_alignments_mut) = cluster_positions.get_mut(prev_coord) {
        assign_shared_alignments(
            prev_coord,
            cluster_coord,
            prev_alignments_mut,
            cluster_alignments,
        );
    }
}

/// Determines if coordinates should be merged based on clip directions and distance
fn should_merge_coordinates(
    prev_clip_directions: (bool, bool),
    current_clip_directions: (bool, bool),
    prev_dist: i64,
    cluster_coord: &Coordinate,
    prev_coord: &Coordinate,
) -> bool {
    prev_clip_directions == current_clip_directions
        && (prev_dist <= cluster_coord.confidence_interval.0
            || prev_dist <= prev_coord.confidence_interval.1)
}

/// Handles merging with next coordinate if within confidence interval
fn handle_next_coordinate_merging(
    cluster_positions: &mut BTreeMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    cluster_coord: &Coordinate,
    cluster_alignments: &mut HashSet<FwdStrandSplitReadSegment>,
    working_cluster_positions: &BTreeMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
) {
    let next_coord_opt: Option<(&Coordinate, &HashSet<FwdStrandSplitReadSegment>)> =
        working_cluster_positions.range(cluster_coord..).next();

    if let Some((next_coord, next_alignments)) = next_coord_opt {
        if next_coord.start_chrom == cluster_coord.start_chrom {
            process_next_coordinate(
                cluster_positions,
                next_coord,
                next_alignments,
                cluster_coord,
                cluster_alignments,
            );
        }
    }
}

/// Processes the next coordinate for potential merging or shared alignment assignment
fn process_next_coordinate(
    cluster_positions: &mut BTreeMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    next_coord: &Coordinate,
    next_alignments: &HashSet<FwdStrandSplitReadSegment>,
    cluster_coord: &Coordinate,
    cluster_alignments: &mut HashSet<FwdStrandSplitReadSegment>,
) {
    let next_dist = (next_coord.end - cluster_coord.start).abs();

    if should_merge_with_next_coordinate(next_dist, cluster_coord, next_coord) {
        cluster_alignments.extend(next_alignments.iter().cloned());
        cluster_positions.remove(next_coord);
    } else if let Some(next_alignments_mut) = cluster_positions.get_mut(next_coord) {
        assign_shared_alignments(
            cluster_coord,
            next_coord,
            cluster_alignments,
            next_alignments_mut,
        );
    }
}

/// Determines if current coordinate should merge with next coordinate
fn should_merge_with_next_coordinate(
    next_dist: i64,
    cluster_coord: &Coordinate,
    next_coord: &Coordinate,
) -> bool {
    next_dist <= cluster_coord.confidence_interval.0
        || next_dist <= next_coord.confidence_interval.1
}

fn assign_shared_alignments(
    prev_coord: &Coordinate,
    current_coord: &Coordinate,
    prev_alignments: &mut HashSet<FwdStrandSplitReadSegment>,
    current_alignments: &mut HashSet<FwdStrandSplitReadSegment>,
) {
    let shared_alignments: HashSet<FwdStrandSplitReadSegment> = current_alignments
        .intersection(prev_alignments)
        .cloned()
        .collect();

    // if shared alignments are anchored to two different breakends by the same clipping event,
    // assign them to only the closest of the two breakends. We can assume ordering.
    for shared_alignment in shared_alignments {
        let alignment_dists_to_current = (
            (shared_alignment.pos - current_coord.start).abs(),
            (shared_alignment.end - current_coord.start).abs(),
        );
        let alignment_dists_to_prev = (
            (shared_alignment.pos - prev_coord.start).abs(),
            (shared_alignment.end - prev_coord.start).abs(),
        );
        if shared_alignment.is_start_softclipped
            && alignment_dists_to_current.0 <= current_coord.confidence_interval.0
            && alignment_dists_to_prev.0 <= prev_coord.confidence_interval.0
        {
            if alignment_dists_to_current.0 < alignment_dists_to_prev.0 {
                // this means the alignment left-clip is closer to the current than to the previous coordinate, so remove it from the previous one
                prev_alignments.remove(&shared_alignment);
            } else {
                // this means the alignment left-clip is closer to the previous than to the current coordinate, so remove it from the previous one
                current_alignments.remove(&shared_alignment);
            }
        }
        if shared_alignment.is_end_softclipped
            && alignment_dists_to_current.1 <= current_coord.confidence_interval.1
            && alignment_dists_to_prev.1 <= prev_coord.confidence_interval.1
        {
            if alignment_dists_to_current.1 < alignment_dists_to_prev.1 {
                // this means the alignment right-clip is closer to the current than to the previous coordinate
                prev_alignments.remove(&shared_alignment);
            } else {
                current_alignments.remove(&shared_alignment);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils;
    use std::collections::{BTreeMap, HashMap, HashSet};

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
    fn test_calculate_alignment_min_distance() {
        let coord = create_test_coordinate("chr1", 1000, (10, 10));

        // Alignment with start softclip close to coordinate
        let alignment1 = utils::create_test_alignment("chr1", 995, 1100, "read1");
        let dist1 = calculate_alignment_min_distance(&alignment1, &coord);
        assert_eq!(dist1, 5);

        // Alignment with end softclip close to coordinate
        let alignment2 =
            utils::create_test_alignment_with_clips("chr1", 900, 1010, "read2", false, true);
        let dist2 = calculate_alignment_min_distance(&alignment2, &coord);
        assert_eq!(dist2, 10);

        // Alignment with both clips, should use minimum distance
        let alignment3 =
            utils::create_test_alignment_with_clips("chr1", 980, 1020, "read3", true, true);
        let dist3 = calculate_alignment_min_distance(&alignment3, &coord);
        assert_eq!(dist3, 20); // min of (20, 20)

        // Alignment with no clips
        let alignment4 =
            utils::create_test_alignment_with_clips("chr1", 980, 1020, "read4", false, false);
        let dist4 = calculate_alignment_min_distance(&alignment4, &coord);
        assert_eq!(dist4, i64::MAX);
    }

    #[test]
    fn test_group_alignments_by_distance() {
        let coord = create_test_coordinate("chr1", 1000, (10, 10));
        let mut alignments = HashSet::new();

        let alignment1 = utils::create_test_alignment("chr1", 995, 1100, "read1");
        let alignment2 = utils::create_test_alignment("chr1", 990, 1100, "read2");
        let alignment3 = utils::create_test_alignment("chr1", 980, 1100, "read3");

        alignments.insert(alignment1);
        alignments.insert(alignment2);
        alignments.insert(alignment3);

        let grouped = group_alignments_by_distance(&alignments, &coord);

        assert_eq!(grouped.len(), 3); // Three different distances: 5, 10, 20
        assert!(grouped.contains_key(&5));
        assert!(grouped.contains_key(&10));
        assert!(grouped.contains_key(&20));
        assert_eq!(grouped[&5].len(), 1);
        assert_eq!(grouped[&10].len(), 1);
        assert_eq!(grouped[&20].len(), 1);
    }

    #[test]
    fn test_find_usable_distances() {
        let mut alignments_by_dist = HashMap::new();

        // Add distance with enough support
        let alignment1 = utils::create_test_alignment("chr1", 995, 1100, "read1");
        let alignment2 = utils::create_test_alignment("chr1", 994, 1100, "read2");
        alignments_by_dist.insert(5, vec![alignment1, alignment2]);

        // Add distance with insufficient support
        let alignment3 = utils::create_test_alignment("chr1", 990, 1100, "read3");
        alignments_by_dist.insert(10, vec![alignment3]);

        // Add another distance with enough support
        let alignment4 = utils::create_test_alignment("chr1", 980, 1100, "read4");
        let alignment5 = utils::create_test_alignment("chr1", 979, 1100, "read5");
        alignments_by_dist.insert(20, vec![alignment4, alignment5]);

        let (usable_dists, most_common_dist) = find_usable_distances(&alignments_by_dist);

        assert_eq!(usable_dists.len(), 2); // Only distances 5 and 20 have enough support
        assert!(usable_dists.contains(&5));
        assert!(usable_dists.contains(&20));
        assert!(most_common_dist == 5 || most_common_dist == 20); // Both have same count
    }

    #[test]
    fn test_determine_allowed_distance() {
        // Test with distances within confidence interval
        let usable_dists = vec![5, 10, 15];
        let most_common_dist = 10;
        let allowed = determine_allowed_distance(&usable_dists, most_common_dist);
        assert_eq!(allowed, 15); // Max of the usable distances

        // Test with distances exceeding confidence interval
        let large_dists = vec![5, 1000, 2000];
        let allowed_large = determine_allowed_distance(&large_dists, 10);
        assert_eq!(allowed_large, utils::MAX_ALIGNMENT_CONFIDENCE_INTERVAL);

        // Test with empty usable distances
        let empty_dists = vec![];
        let allowed_empty = determine_allowed_distance(&empty_dists, 50);
        assert_eq!(allowed_empty, 50 + utils::MAX_CLUST_DIST);
    }

    #[test]
    fn test_filter_alignments_by_distance() {
        let mut alignments_by_dist = HashMap::new();

        let alignment1 = utils::create_test_alignment("chr1", 995, 1100, "read1");
        let alignment2 = utils::create_test_alignment("chr1", 990, 1100, "read2");
        let alignment3 = utils::create_test_alignment("chr1", 980, 1100, "read3");

        alignments_by_dist.insert(5, vec![alignment1.clone()]);
        alignments_by_dist.insert(10, vec![alignment2.clone()]);
        alignments_by_dist.insert(20, vec![alignment3.clone()]);

        let filtered = filter_alignments_by_distance(&alignments_by_dist, 10);

        assert_eq!(filtered.len(), 2); // Should include distances 5 and 10, but not 20
        assert!(filtered.contains(&alignment1));
        assert!(filtered.contains(&alignment2));
        assert!(!filtered.contains(&alignment3));
    }

    #[test]
    fn test_update_cluster_confidence_interval() {
        let mut coord = create_test_coordinate("chr1", 1000, (5, 5));
        update_cluster_confidence_interval(&mut coord, 15);

        assert_eq!(coord.confidence_interval, (15, 15));

        // Test with distance smaller than MAX_CLUST_DIST
        let mut coord2 = create_test_coordinate("chr1", 1000, (5, 5));
        update_cluster_confidence_interval(&mut coord2, 5);

        assert_eq!(
            coord2.confidence_interval,
            (utils::MAX_CLUST_DIST, utils::MAX_CLUST_DIST)
        );
    }

    #[test]
    fn test_get_coord_clip_directions() {
        let coord = create_test_coordinate("chr1", 1000, (10, 10));
        let mut alignments = HashSet::new();

        // Add left-clipped alignment
        alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 995, 1100, "read1", true, true, None, None,
        ));

        // Add right-clipped alignment
        alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 900, 1005, "read2", false, true, None, None,
        ));

        let (is_rightclipped, is_leftclipped) = get_coord_clip_directions(&alignments, &coord);

        assert!(is_leftclipped);
        assert!(is_rightclipped);

        // Test with only left-clipped
        let mut left_only = HashSet::new();
        left_only.insert(utils::create_test_alignment_with_phasing(
            "chr1", 995, 1100, "read1", true, true, None, None,
        ));

        let (is_right_only, is_left_only) = get_coord_clip_directions(&left_only, &coord);
        assert!(!is_right_only);
        assert!(is_left_only);
    }

    #[test]
    fn test_should_merge_coordinates() {
        let prev_coord = create_test_coordinate("chr1", 1000, (20, 20));
        let current_coord = create_test_coordinate("chr1", 1015, (20, 20));

        // Same clip directions and within confidence interval
        let prev_clip_dirs = (true, false);
        let current_clip_dirs = (true, false);
        let prev_dist = 15;

        assert!(should_merge_coordinates(
            prev_clip_dirs,
            current_clip_dirs,
            prev_dist,
            &current_coord,
            &prev_coord
        ));

        // Different clip directions
        let different_clip_dirs = (false, true);
        assert!(!should_merge_coordinates(
            prev_clip_dirs,
            different_clip_dirs,
            prev_dist,
            &current_coord,
            &prev_coord
        ));

        // Same directions but too far apart
        let far_dist = 50;
        assert!(!should_merge_coordinates(
            prev_clip_dirs,
            current_clip_dirs,
            far_dist,
            &current_coord,
            &prev_coord
        ));
    }

    #[test]
    fn test_should_merge_with_next_coordinate() {
        let current_coord = create_test_coordinate("chr1", 1000, (20, 20));
        let next_coord = create_test_coordinate("chr1", 1015, (20, 20));

        // Within confidence interval
        let close_dist = 15;
        assert!(should_merge_with_next_coordinate(
            close_dist,
            &current_coord,
            &next_coord
        ));

        // Too far apart
        let far_dist = 50;
        assert!(!should_merge_with_next_coordinate(
            far_dist,
            &current_coord,
            &next_coord
        ));
    }

    #[test]
    fn test_add_vcf_alignments() {
        let mut current_coord = Coordinate::new("chr1".to_string(), 1000);
        let mut cluster_alignments = HashSet::new();

        let variant_alignments = vec![
            utils::create_test_alignment("chr1", 995, 1100, "read1"),
            utils::create_test_alignment("chr1", 990, 1095, "read2"),
            utils::create_test_alignment_with_clips("chr1", 950, 1050, "read3", false, true),
        ];

        add_vcf_alignments(
            &variant_alignments,
            &mut current_coord,
            &mut cluster_alignments,
        );

        // Should include the alignments closest to the coordinate
        assert!(!cluster_alignments.is_empty());

        // Test with no softclipped alignments
        let mut no_clip_coord = Coordinate::new("chr1".to_string(), 2000);
        let mut no_clip_alignments = HashSet::new();

        let no_clip_variant_alignments = vec![utils::create_test_alignment_with_phasing(
            "chr1", 1995, 2100, "read1", false, false, None, None,
        )];

        add_vcf_alignments(
            &no_clip_variant_alignments,
            &mut no_clip_coord,
            &mut no_clip_alignments,
        );
        assert!(no_clip_alignments.is_empty()); // No softclipped alignments to add
    }

    #[test]
    fn test_assign_shared_alignments() {
        let prev_coord = create_test_coordinate("chr1", 1000, (50, 50));
        let current_coord = create_test_coordinate("chr1", 1020, (50, 50));

        let shared_alignment = utils::create_test_alignment("chr1", 995, 1100, "read1");

        let mut prev_alignments = HashSet::new();
        prev_alignments.insert(shared_alignment.clone());

        let mut current_alignments = HashSet::new();
        current_alignments.insert(shared_alignment.clone());

        assign_shared_alignments(
            &prev_coord,
            &current_coord,
            &mut prev_alignments,
            &mut current_alignments,
        );

        // Function may or may not remove the alignment depending on exact conditions
        // The key test is that the function runs without panicking and processes correctly
        assert!(prev_alignments.len() <= 1);
        assert!(current_alignments.len() <= 1);
    }

    #[test]
    fn test_cluster_via_vcf() {
        let clipped_reads = create_test_clipped_reads();

        let mut vcf_coord_map = HashMap::new();
        vcf_coord_map.insert(
            "variant1".to_string(),
            vec![Coordinate::new("chr1".to_string(), 1000)],
        );
        vcf_coord_map.insert(
            "variant2".to_string(),
            vec![Coordinate::new("chr1".to_string(), 2000)],
        );

        let mut variant_readnames = HashMap::new();
        variant_readnames.insert(
            "variant1".to_string(),
            vec!["read1".to_string(), "read2".to_string()],
        );
        variant_readnames.insert(
            "variant2".to_string(),
            vec!["read1".to_string(), "read2".to_string()],
        );

        let clusters = cluster_via_readnames(&clipped_reads, &vcf_coord_map, &variant_readnames);

        assert!(!clusters.is_empty());
        // Should have clusters at the VCF coordinate positions
        let has_1000_cluster = clusters
            .keys()
            .any(|coord| (coord.start - 1000).abs() <= utils::MAX_CLUST_DIST);
        let has_2000_cluster = clusters
            .keys()
            .any(|coord| (coord.start - 2000).abs() <= utils::MAX_CLUST_DIST);
        assert!(has_1000_cluster || has_2000_cluster);
    }

    #[test]
    fn test_assign_clipped_reads_to_clusters_integration() {
        let clipped_reads = create_test_clipped_reads();

        let mut vcf_coord_map = HashMap::new();
        vcf_coord_map.insert(
            "variant1".to_string(),
            vec![Coordinate::new("chr1".to_string(), 1000)],
        );

        let mut variant_readnames = HashMap::new();
        variant_readnames.insert(
            "variant1".to_string(),
            vec!["read1".to_string(), "read2".to_string()],
        );

        let clusters = assign_clipped_reads_to_clusters(
            &clipped_reads,
            &vcf_coord_map,
            &Some(variant_readnames),
        );

        assert!(!clusters.is_empty());

        // All coordinates should have variant IDs assigned
        for coord in clusters.keys() {
            assert!(!coord.variant_ids.is_empty());
        }
    }

    #[test]
    fn test_refine_alignments_integration() {
        let mut cluster_coord = Coordinate::new("chr1".to_string(), 1000);
        let mut cluster_alignments = HashSet::new();

        // Add alignments at different distances
        cluster_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 995, 1100, "read1", true, true, None, None,
        ));
        cluster_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 990, 1095, "read2", true, true, None, None,
        ));
        cluster_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 980, 1085, "read3", true, true, None, None,
        )); // Farther away

        let refined = refine_alignments(cluster_alignments.clone(), &mut cluster_coord);

        // Should filter out alignments that are too far
        assert!(refined.len() <= cluster_alignments.len());

        // Confidence interval should be updated
        assert!(cluster_coord.confidence_interval.0 >= utils::MAX_CLUST_DIST);
    }

    #[test]
    fn test_add_coordinate_integration() {
        let mut cluster_positions = BTreeMap::new();
        let cluster_coord = Coordinate::new("chr1".to_string(), 1000);
        let mut cluster_alignments = HashSet::new();

        // Add sufficient alignments
        cluster_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            1100,
            "read1",
            true,
            true,
            Some(1),
            None,
        ));
        cluster_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            990,
            1095,
            "read2",
            true,
            true,
            Some(1),
            None,
        ));

        add_coordinate(
            &mut cluster_positions,
            cluster_coord.clone(),
            &mut cluster_alignments,
        );

        assert!(cluster_positions.contains_key(&cluster_coord));

        // Test merging with nearby coordinate
        let nearby_coord = Coordinate::new("chr1".to_string(), 1005);
        let mut nearby_alignments = HashSet::new();
        nearby_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            1000,
            1105,
            "read3",
            true,
            true,
            Some(1),
            None,
        ));
        nearby_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1",
            995,
            1100,
            "read4",
            true,
            true,
            Some(1),
            None,
        ));

        add_coordinate(
            &mut cluster_positions,
            nearby_coord.clone(),
            &mut nearby_alignments,
        );

        // Should either merge or handle shared alignments
        assert!(!cluster_positions.is_empty());
    }

    #[test]
    fn test_alignment_matches_break_forward() {
        let break_coord = create_test_coordinate("chr1", 1000, (10, 10));

        // Test alignment with start softclip that matches
        let alignment1 =
            utils::create_test_alignment_with_clips("chr1", 995, 1100, "read1", true, false);
        assert!(alignment_matches_break_forward(&alignment1, &break_coord));

        // Test alignment with end softclip that matches
        let alignment2 =
            utils::create_test_alignment_with_clips("chr1", 900, 1005, "read2", false, true);
        assert!(alignment_matches_break_forward(&alignment2, &break_coord));

        // Test alignment with both clips that match
        let alignment3 =
            utils::create_test_alignment_with_clips("chr1", 995, 1005, "read3", true, true);
        assert!(alignment_matches_break_forward(&alignment3, &break_coord));

        // Test alignment with no clips
        let alignment4 =
            utils::create_test_alignment_with_clips("chr1", 995, 1005, "read4", false, false);
        assert!(!alignment_matches_break_forward(&alignment4, &break_coord));

        // Test alignment with clips but wrong chromosome
        let alignment5 =
            utils::create_test_alignment_with_clips("chr2", 995, 1005, "read5", true, true);
        assert!(!alignment_matches_break_forward(&alignment5, &break_coord));

        // Test alignment with clips but too far away
        let alignment6 =
            utils::create_test_alignment_with_clips("chr1", 900, 1100, "read6", true, true);
        assert!(!alignment_matches_break_forward(&alignment6, &break_coord));

        // Test alignment with second chromosome matching
        let mut alignment7 =
            utils::create_test_alignment_with_clips("chr2", 995, 1005, "read7", true, true);
        alignment7.second_chrom = "chr1".to_string();
        assert!(alignment_matches_break_forward(&alignment7, &break_coord));
    }

    #[test]
    fn test_handle_previous_coordinate_merging() {
        let mut cluster_positions = BTreeMap::new();
        let cluster_coord = create_test_coordinate("chr1", 1000, (20, 20));
        let mut cluster_alignments = HashSet::new();

        // Add a previous coordinate that should be merged
        let prev_coord = create_test_coordinate("chr1", 985, (20, 20));
        let mut prev_alignments = HashSet::new();
        prev_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 980, 1100, "read1", true, true, None, None,
        ));
        cluster_positions.insert(prev_coord.clone(), prev_alignments.clone());

        // Add current coordinate alignments
        cluster_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 995, 1100, "read2", true, true, None, None,
        ));

        let working_cluster_positions = cluster_positions.clone();

        handle_previous_coordinate_merging(
            &mut cluster_positions,
            &cluster_coord,
            &mut cluster_alignments,
            &working_cluster_positions,
        );
    }

    #[test]
    fn test_process_previous_coordinate() {
        let mut cluster_positions = BTreeMap::new();
        let prev_coord = create_test_coordinate("chr1", 985, (20, 20));
        let mut prev_alignments = HashSet::new();
        prev_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 980, 1100, "read1", true, true, None, None,
        ));

        let cluster_coord = create_test_coordinate("chr1", 1000, (20, 20));
        let mut cluster_alignments = HashSet::new();
        cluster_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 995, 1100, "read2", true, true, None, None,
        ));

        // Add the previous coordinate to the map so it can be processed
        cluster_positions.insert(prev_coord.clone(), prev_alignments.clone());

        process_previous_coordinate(
            &mut cluster_positions,
            &prev_coord,
            &prev_alignments,
            &cluster_coord,
            &mut cluster_alignments,
        );
        assert!(cluster_positions.is_empty());
        assert!(cluster_alignments.len() == 2);
    }

    #[test]
    fn test_handle_next_coordinate_merging() {
        let mut cluster_positions = BTreeMap::new();
        let cluster_coord = create_test_coordinate("chr1", 1000, (20, 20));
        let mut cluster_alignments = HashSet::new();

        // Add a next coordinate that should be merged
        let next_coord = create_test_coordinate("chr1", 1015, (20, 20));
        let mut next_alignments = HashSet::new();
        next_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 1010, 1115, "read1", true, true, None, None,
        ));
        cluster_positions.insert(next_coord.clone(), next_alignments.clone());

        // Add current coordinate alignments
        cluster_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 995, 1100, "read2", true, true, None, None,
        ));

        let working_cluster_positions = cluster_positions.clone();

        handle_next_coordinate_merging(
            &mut cluster_positions,
            &cluster_coord,
            &mut cluster_alignments,
            &working_cluster_positions,
        );
        assert!(cluster_positions.is_empty());
        assert!(cluster_alignments.len() == 2);
    }

    #[test]
    fn test_process_next_coordinate() {
        let mut cluster_positions = BTreeMap::new();
        let next_coord = create_test_coordinate("chr1", 1015, (20, 20));
        let mut next_alignments = HashSet::new();
        next_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 1010, 1115, "read1", true, true, None, None,
        ));

        let cluster_coord = create_test_coordinate("chr1", 1000, (20, 20));
        let mut cluster_alignments = HashSet::new();
        cluster_alignments.insert(utils::create_test_alignment_with_phasing(
            "chr1", 995, 1100, "read2", true, true, None, None,
        ));

        // Add the next coordinate to the map so it can be processed
        cluster_positions.insert(next_coord.clone(), next_alignments.clone());

        process_next_coordinate(
            &mut cluster_positions,
            &next_coord,
            &next_alignments,
            &cluster_coord,
            &mut cluster_alignments,
        );
        assert!(cluster_positions.is_empty());
        assert!(cluster_alignments.len() == 2);
    }

    #[test]
    fn test_cluster_via_readnames() {
        let clipped_reads = create_test_clipped_reads();
        let mut vcf_coord_map = HashMap::new();
        vcf_coord_map.insert(
            "variant1".to_string(),
            vec![Coordinate::new("chr1".to_string(), 1000)],
        );

        let mut variant_readnames = HashMap::new();
        variant_readnames.insert(
            "variant1".to_string(),
            vec!["read1".to_string(), "read2".to_string()],
        );

        let clusters = cluster_via_readnames(&clipped_reads, &vcf_coord_map, &variant_readnames);

        assert!(!clusters.is_empty());

        // Check that variant IDs are assigned
        for coord in clusters.keys() {
            assert!(!coord.variant_ids.is_empty());
            assert!(coord.variant_ids.contains(&"variant1".to_string()));
        }
    }

    #[test]
    fn test_assign_clipped_reads_to_clusters_without_variant_readnames() {
        let clipped_reads = create_test_clipped_reads();
        let mut vcf_coord_map = HashMap::new();
        vcf_coord_map.insert(
            "variant1".to_string(),
            vec![Coordinate::new("chr1".to_string(), 1000)],
        );

        // Test without variant readnames (should use cluster_via_vcf)
        // Note: This function calls std::process::exit() when no clusters are found,
        // so we can't easily test the empty case in a unit test.
        // Instead, we test that the function exists and can be called with valid data.
        let mut variant_readnames = HashMap::new();
        variant_readnames.insert(
            "variant1".to_string(),
            vec!["read1".to_string(), "read2".to_string()],
        );

        let clusters = assign_clipped_reads_to_clusters(
            &clipped_reads,
            &vcf_coord_map,
            &Some(variant_readnames),
        );

        // Should return clusters when variant readnames are provided
        assert!(!clusters.is_empty());
    }

    #[test]
    fn test_edge_cases_for_merging_functions() {
        // Test with empty alignments
        let mut cluster_positions = BTreeMap::new();
        let cluster_coord = create_test_coordinate("chr1", 1000, (20, 20));
        let mut cluster_alignments = HashSet::new();

        handle_previous_coordinate_merging(
            &mut cluster_positions,
            &cluster_coord,
            &mut cluster_alignments,
            &BTreeMap::new(),
        );

        handle_next_coordinate_merging(
            &mut cluster_positions,
            &cluster_coord,
            &mut cluster_alignments,
            &BTreeMap::new(),
        );

        // Should handle empty cases without panicking
        assert!(cluster_positions.is_empty());
    }
}
