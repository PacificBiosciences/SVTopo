use log::{debug, error};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::time::SystemTime;
use std::{i64, u32};

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

/// Given a hashmap with groups of read alignments keyed by readnames,
/// identify all clusters of softclip locations with at least two supporting reads.
/// These will define the start and stop coordinate of blocks in complex SVs.
/// Join together the softclips that are within max_clust_width bases of one another and
/// collect these in a hashmap of coordinates to vectors of supporting alignments.
pub fn assign_clipped_reads_to_clusters(
    clipped_reads: &HashMap<String, Vec<FwdStrandSplitReadSegment>>,
    vcf_coord_map: &HashMap<String, Vec<Coordinate>>,
    variant_readnames: &HashMap<String, Vec<String>>,
) -> HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>> {
    let start_time = SystemTime::now();
    let cluster_positions = cluster_via_vcf(clipped_reads, vcf_coord_map, variant_readnames);

    if cluster_positions.is_empty() {
        let msg1 = "No readnames were matched to variant IDs from the VCF.";
        let msg2 = "This could be caused several things, including a file mismatch such as between the BAM and the VCF, or the VCF and JSON.";

        error!("{} {}", msg1, msg2);
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
fn cluster_via_vcf(
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

/// Add alignments to the cluster alignments, from vcf coordinate matching
fn add_vcf_alignments(
    variant_alignments: &[FwdStrandSplitReadSegment],
    current_coord: &mut Coordinate,
    cluster_alignments: &mut HashSet<FwdStrandSplitReadSegment>,
) {
    let mut closest_alignment_dist = i64::max_value();

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
    let mut refined_cluster_alignments: HashSet<FwdStrandSplitReadSegment> = HashSet::new();
    if cluster_alignments.is_empty() {
        return cluster_alignments;
    }
    let mut alignments_by_dist = HashMap::new();
    let mut ordered_cluster_alignments: Vec<&FwdStrandSplitReadSegment> =
        cluster_alignments.iter().collect();
    ordered_cluster_alignments.sort();
    for alignment in ordered_cluster_alignments {
        let mut left_dist = std::i64::MAX;
        let mut right_dist = std::i64::MAX;
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
        alignments_by_dist
            .entry(std::cmp::min(left_dist, right_dist))
            .or_insert_with(Vec::new)
            .push(alignment.clone());
    }
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

    // if the range of clipping locations that support this breakend is larger than allowed, limit
    // to only the aligments that are closer than that distance.
    let mut allowed_alignment_distance_opt = None;
    if let (Some(leftmost), Some(rightmost)) = (
        usable_clip_site_dists.first(),
        usable_clip_site_dists.last(),
    ) {
        if rightmost - leftmost > utils::MAX_ALIGNMENT_CONFIDENCE_INTERVAL {
            allowed_alignment_distance_opt = Some(utils::MAX_ALIGNMENT_CONFIDENCE_INTERVAL);
        } else {
            allowed_alignment_distance_opt = Some(std::cmp::max(*rightmost, *leftmost));
        }
    }
    if allowed_alignment_distance_opt.is_none() {
        allowed_alignment_distance_opt = Some(most_common_dist);
        if let Some(dist) = allowed_alignment_distance_opt {
            allowed_alignment_distance_opt = Some(dist + utils::MAX_CLUST_DIST);
        }
    }
    if let Some(allowed_alignment_distance) = allowed_alignment_distance_opt {
        for (dist, alignments) in alignments_by_dist {
            if dist <= allowed_alignment_distance {
                refined_cluster_alignments.extend(alignments.into_iter());
            }
        }
        let ci = std::cmp::max(utils::MAX_CLUST_DIST, allowed_alignment_distance);
        cluster_coord.confidence_interval = (ci, ci);
    }
    refined_cluster_alignments
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
    let working_cluster_positions = cluster_positions.clone();
    if cluster_alignments.len() >= utils::MIN_CONNECTION_SUPPORT as usize {
        let prev_coord_opt: Option<(&Coordinate, &HashSet<FwdStrandSplitReadSegment>)> =
            working_cluster_positions
                .range(..&cluster_coord)
                .next_back();
        if let Some((prev_coord, prev_alignments)) = prev_coord_opt {
            if prev_coord.end_chrom == cluster_coord.start_chrom {
                let prev_clip_directions = get_coord_clip_directions(prev_alignments, prev_coord);
                let current_clip_directions =
                    get_coord_clip_directions(cluster_alignments, &cluster_coord);
                let prev_dist = (prev_coord.end - cluster_coord.start).abs();

                if prev_clip_directions == current_clip_directions
                    && (prev_dist <= cluster_coord.confidence_interval.0
                        || prev_dist <= prev_coord.confidence_interval.1)
                {
                    cluster_alignments.extend(prev_alignments.iter().cloned());
                    cluster_positions.remove_entry(prev_coord);
                } else if let Some(prev_alignments_mut) = cluster_positions.get_mut(prev_coord) {
                    // make sure that aligment ends are only associated with the best-match coordinate
                    assign_shared_alignments(
                        prev_coord,
                        &cluster_coord,
                        prev_alignments_mut,
                        cluster_alignments,
                    );
                }
            }
        }
        let next_coord_opt: Option<(&Coordinate, &HashSet<FwdStrandSplitReadSegment>)> =
            working_cluster_positions.range(&cluster_coord..).next();
        if let Some((next_coord, next_alignments)) = next_coord_opt {
            if next_coord.start_chrom == cluster_coord.start_chrom {
                let next_dist = (next_coord.end - cluster_coord.start).abs();
                if next_dist <= cluster_coord.confidence_interval.0
                    || next_dist <= next_coord.confidence_interval.1
                {
                    cluster_alignments.extend(next_alignments.iter().cloned());
                    cluster_positions.remove(next_coord);
                } else if let Some(next_alignments_mut) = cluster_positions.get_mut(next_coord) {
                    // make sure that aligment ends are only associated with the best-match coordinate
                    assign_shared_alignments(
                        &cluster_coord,
                        next_coord,
                        cluster_alignments,
                        next_alignments_mut,
                    );
                }
            }
        }
        cluster_positions.insert(cluster_coord, cluster_alignments.clone());
    }
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

/// Combine entries of softclip locations where multiple softclips are very close together
/// by storing a position, walking through the list linearly, and storing the cluster when
/// a new position is processed that's more than max_clust_width away
fn cluster_via_clipping(
    clipped_reads: &HashMap<String, Vec<FwdStrandSplitReadSegment>>,
    allow_unphased: bool,
) -> HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>> {
    let mut clip_coords = HashMap::new();

    // add all softclip locations to a hashmap of coordinate to name of supporting read
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
    debug!(
        "{} putative break locations pre-filtering",
        clip_coords.len()
    );
    for (break_location, support) in clip_coords.iter() {
        debug!("{} (support={})", break_location, support.len());
    }
    let mut clip_coords_vec: Vec<Coordinate> = clip_coords.keys().cloned().collect();
    clip_coords_vec.sort();
    let mut cluster_positions = HashMap::new();
    let mut begin_iter: usize = 0;
    for (end_iter, coord) in clip_coords_vec.iter().enumerate() {
        let cluster_start = &clip_coords_vec[begin_iter];
        // if the current coordinate is too far from the start of the cluster,
        // store the cluster and move on
        if (coord.start_chrom != cluster_start.start_chrom)
            || (coord.start > (cluster_start.start + cluster_start.confidence_interval.0))
        {
            let middle = get_clust_middle(&clip_coords_vec[begin_iter..end_iter]);
            let alignments =
                get_clust_aligments(&clip_coords_vec[begin_iter..end_iter], &clip_coords);
            let spanned_count = alignments.iter().filter(|a| a.spans).count();
            if spanned_count >= utils::MIN_CONNECTION_SUPPORT as usize {
                add_cluster(
                    middle,
                    alignments,
                    &mut cluster_positions,
                    utils::MIN_CONNECTION_SUPPORT,
                    allow_unphased,
                );
            }
            begin_iter = end_iter;
        }
    }
    if !cluster_positions.is_empty() {
        //add the last cluster
        let middle = get_clust_middle(&clip_coords_vec[begin_iter..]);
        let alignments = get_clust_aligments(&clip_coords_vec[begin_iter..], &clip_coords);
        add_cluster(
            middle,
            alignments,
            &mut cluster_positions,
            utils::MIN_CONNECTION_SUPPORT,
            allow_unphased,
        );
    }
    cluster_positions
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
