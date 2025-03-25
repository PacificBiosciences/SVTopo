use crate::{
    containers::{ComplexSVCalls, CoverageMap},
    utils,
};
use log::debug;
use std::{path::PathBuf, time::SystemTime};

/// Loop through the BAM for all genomic intervals that have been connected.
/// If the coverage is within expected thresholds and the mapq is good,
/// keep the blocks. Else, drop them from return.
pub fn apply_filters(
    annotated_graphs: &ComplexSVCalls,
    bam_filename: PathBuf,
    max_allowed_coverage: u32,
) -> ComplexSVCalls {
    let start_time = SystemTime::now();
    let mut coverage_map: CoverageMap = CoverageMap::new(&annotated_graphs.event_graphs);
    debug!("Created coverage map");
    coverage_map.update_coverages(bam_filename);
    debug!("Update coverage map from bam");

    let mut filtered_annotated_graphs = Vec::new();
    for annotated_graph in annotated_graphs.event_graphs.iter() {
        let mut good_coverage = true;
        for block in annotated_graph.iter() {
            debug!("{}", block.region);
            if let Some((covs, low_mapq_covs)) = coverage_map.coverages.get(&block.region) {
                if let Some(max_cov) = covs.iter().max() {
                    debug!("{} max coverage: {}", block.region, max_cov,);
                    if (*max_cov) > max_allowed_coverage {
                        good_coverage = false;
                        debug!(
                            "{} skipped: Max coverage too high ({}), above {}",
                            block.region, max_cov, max_allowed_coverage,
                        );
                    }
                }
                let cov_sum: u32 = covs.iter().sum();
                if cov_sum > 0 {
                    let low_mapq_cov_sum: u32 = low_mapq_covs.iter().sum();
                    let low_mapq_fract = (low_mapq_cov_sum / cov_sum) as f64;
                    debug!(
                        "{} low MAPQ coverage fraction: {}",
                        block.region, low_mapq_fract
                    );
                    if low_mapq_fract > utils::EXTREMELY_LOW_MAPQ_FRACTION_THRESHOLD {
                        good_coverage = false;
                        debug!(
                            "{} skipped: Low MAPQ fraction too high ({}), above {}",
                            block.region,
                            low_mapq_fract,
                            utils::EXTREMELY_LOW_MAPQ_FRACTION_THRESHOLD
                        );
                    }
                }
            }
        }
        if good_coverage {
            filtered_annotated_graphs.push(annotated_graph.clone());
        }
    }

    debug!(
        "Filtering final events: {}s",
        start_time.elapsed().unwrap().as_secs()
    );
    ComplexSVCalls::new(filtered_annotated_graphs)
}
