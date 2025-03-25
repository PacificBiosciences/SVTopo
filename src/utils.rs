use std::io::{BufRead, BufReader};

use log::error;

/// Greatest distance in bp allowed between softclipping
/// locations to treat them as the same clipping event.
pub const MAX_CLUST_DIST: i64 = 15;

/// Minumum number of alignment connections required to
/// treat a pair of softclipping locations as connected
pub const MIN_CONNECTION_SUPPORT: u32 = 2;

/// Maximum distance on the same chromosome to join via phasing
pub const MAX_PHASE_DIST: i64 = 500_000;

/// Maximum number of blocks allowed in a single event graph
pub const MAX_GRAPH_SIZE: usize = 20;

/// Minimum amount of clipping to consider this read as a clipped
/// read, supporting a structural rearrangement
pub const MIN_CLIPPING: i64 = 100;

/// Minimum required mapping quality (MAPQ) for a read to be used
pub const MIN_MAPQ: u8 = 20;

/// Allowable fraction of extremely low mapq coverage
pub const EXTREMELY_LOW_MAPQ_FRACTION_THRESHOLD: f64 = 0.05;

/// Definition of extremely low mapq
pub const EXTREMELY_LOW_MAPQ: u8 = 5;

/// Maximum allowed read alignment confidence interval to match sawfish call
pub const MAX_ALIGNMENT_CONFIDENCE_INTERVAL: i64 = 100;

/// first two bytes of a gzip file that indicatee the compression algorithm used
const GZIP_INDICATOR: [u8; 2] = [0x1F, 0x8B];

/// Check a coordinate entry against the exclude regions to see if it should be skipped.
pub fn entry_excluded(
    exclude_regions: &crate::containers::ExcludeRegions,
    entry: (String, i64, i64),
    target_region: &Option<crate::containers::TargetCoordinate>,
) -> bool {
    let mut excluded = false;
    if let Some(exclude_coords) = exclude_regions.regions.get(&entry.0) {
        let start_coord = crate::containers::Coordinate::new(entry.0.clone(), entry.1);
        let end_coord = crate::containers::Coordinate::new(entry.0.clone(), entry.2);
        excluded = exclude_coords.range(start_coord..end_coord).count() > 0;
        if !excluded {
            if let Some(target_coordinate) = target_region {
                let wrong_chrom = entry.0 != target_coordinate.coord.start_chrom;
                let before_target = entry.1 < target_coordinate.coord.start
                    && entry.2 < target_coordinate.coord.start;
                let after_target =
                    entry.1 > target_coordinate.coord.end && entry.2 > target_coordinate.coord.end;
                if wrong_chrom || before_target || after_target {
                    excluded = true;
                }
            }
        }
    }
    excluded
}

pub fn is_local_file(filepath: &String) -> bool {
    let path = std::path::Path::new(filepath);

    match std::fs::metadata(path) {
        Ok(metadata) => metadata.is_file(),
        Err(_) => false, // If there is an error (e.g., path doesn't exist), return false
    }
}

/// Check if a file is a gzipped file from a String path
pub fn is_gzipped(path: &String) -> bool {
    if !is_local_file(path) {
        return false;
    }
    let file_handle = match std::fs::File::open(path) {
        Ok(file) => file,
        Err(e) => {
            error!("File does not exist: \"{}\"", path);
            error!("{}", e);
            std::process::exit(exitcode::IOERR);
        }
    };
    let mut reader = std::io::BufReader::new(file_handle);
    let mut gzip_indicator_bytes = [0; 2];
    let _ = std::io::Read::read_exact(&mut reader, &mut gzip_indicator_bytes);
    let _ = std::io::Seek::rewind(&mut reader);
    gzip_indicator_bytes == GZIP_INDICATOR
}

/// Read a plain text or gzipped text file into vector of Strings by line
pub fn read_file_from_path(file_path: &String) -> Vec<String> {
    assert!(is_local_file(file_path));

    let path = std::path::Path::new(file_path);
    let file: std::fs::File = match std::fs::File::open(path) {
        Ok(file) => file,
        Err(_) => {
            error!("File not found {}", file_path);
            std::process::exit(exitcode::IOERR);
        }
    };

    let lines_result = match is_gzipped(file_path) {
        true => {
            let bgzf_reader = match rust_htslib::bgzf::Reader::from_path(file_path) {
                Ok(r) => r,
                Err(_) => {
                    error!("Failed to read exclude regions file {}", file_path);
                    std::process::exit(exitcode::IOERR);
                }
            };
            let reader = BufReader::new(bgzf_reader);
            reader.lines().collect()
        }
        false => {
            let reader = BufReader::new(file);
            reader.lines().collect()
        }
    };
    match lines_result {
        Ok(l) => l,
        Err(_) => {
            error!("Failed to read exclude regions file {}", file_path);
            std::process::exit(exitcode::IOERR);
        }
    }
}
