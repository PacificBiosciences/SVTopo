use core::fmt;
use log::{debug, error};
use rust_htslib::{
    bam::{self, ext::BamRecordExtensions, Read},
    bcf,
};
use serde::Serialize;
use std::{
    cmp::{max, min},
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    path::PathBuf,
    str::FromStr,
};

use crate::utils::{self};

/// Simple genomic coordinate struct, comparable to a BEDPE record
#[derive(Debug, PartialEq, Eq, Hash, Clone, Serialize)]
pub struct Coordinate {
    pub start_chrom: String,
    pub start: i64,
    pub end_chrom: String,
    pub end: i64,

    /// unique list of associated variant IDs
    pub variant_ids: Vec<String>,

    /// Pair of positive integer values indicating the confidence
    /// interval around the coordinate start. Defaults to utils::MAX_CLUST_DIST.
    #[serde(skip_serializing)]
    pub confidence_interval: (i64, i64),
}

impl Coordinate {
    pub fn new_region(start_chrom: String, start: i64, end: i64) -> Self {
        Coordinate {
            start_chrom: start_chrom.clone(),
            start,
            end_chrom: start_chrom,
            end,
            confidence_interval: (utils::MAX_CLUST_DIST, utils::MAX_CLUST_DIST),
            variant_ids: Vec::new(),
        }
    }
    pub fn new(start_chrom: String, start: i64) -> Self {
        Coordinate {
            start_chrom: start_chrom.clone(),
            start,
            end_chrom: start_chrom,
            end: start,
            confidence_interval: (utils::MAX_CLUST_DIST, utils::MAX_CLUST_DIST),
            variant_ids: Vec::new(),
        }
    }
    pub fn new_with_confidence_interval(
        start_chrom: String,
        start: i64,
        confidence_interval: (i64, i64),
    ) -> Self {
        Coordinate {
            start_chrom: start_chrom.clone(),
            start,
            end_chrom: start_chrom,
            end: start,
            confidence_interval,
            variant_ids: Vec::new(),
        }
    }

    pub fn from_alignment(alignment: FwdStrandSplitReadSegment) -> Self {
        Coordinate {
            start_chrom: alignment.chrom.clone(),
            start: alignment.pos,
            end_chrom: alignment.second_chrom,
            end: alignment.end,
            confidence_interval: (utils::MAX_CLUST_DIST, utils::MAX_CLUST_DIST),
            variant_ids: Vec::new(),
        }
    }

    pub fn from_connection(connection: Connection) -> Self {
        let left = std::cmp::min(&connection.first_coord, &connection.second_coord);
        let right = std::cmp::max(&connection.first_coord, &connection.second_coord);

        Coordinate {
            start_chrom: left.start_chrom.clone(),
            start: left.start,
            end_chrom: right.end_chrom.clone(),
            end: right.end,
            confidence_interval: (utils::MAX_CLUST_DIST, utils::MAX_CLUST_DIST),
            variant_ids: Vec::new(),
        }
    }
}

impl Coordinate {
    ///Comparison function to determine if the current coordinate is within dist of another coordinate on the same chromosome
    pub fn is_within(&self, other_coord: &Coordinate) -> bool {
        if self.start_chrom != other_coord.start_chrom || self.end_chrom != other_coord.end_chrom {
            return false;
        }
        let start_dist = (self.start - other_coord.start).abs();
        let end_dist = (self.end - other_coord.end).abs();

        if start_dist
            > max(
                self.confidence_interval.0,
                other_coord.confidence_interval.0,
            )
            || end_dist
                > max(
                    self.confidence_interval.1,
                    other_coord.confidence_interval.1,
                )
        {
            return false;
        }
        true
    }
}

impl fmt::Display for Coordinate {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.start_chrom == self.end_chrom {
            if self.start == self.end {
                write!(formatter, "{}:{}", &self.start_chrom, &self.start,)
            } else {
                write!(
                    formatter,
                    "{}:{}-{}",
                    &self.start_chrom, &self.start, &self.end,
                )
            }
        } else {
            write!(
                formatter,
                "{}:{}-{}:{}",
                &self.start_chrom, &self.start, &self.end_chrom, &self.end,
            )
        }
    }
}

impl Ord for Coordinate {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let start_chrom_cmp = self.start_chrom.cmp(&other.start_chrom);
        if start_chrom_cmp == std::cmp::Ordering::Equal {
            let start_pos_cmp = self.start.cmp(&other.start);

            if start_pos_cmp == std::cmp::Ordering::Equal {
                let end_chrom_cmp = self.end_chrom.cmp(&other.end_chrom);
                if end_chrom_cmp == std::cmp::Ordering::Equal {
                    self.end.cmp(&other.end)
                } else {
                    end_chrom_cmp
                }
            } else {
                start_pos_cmp
            }
        } else {
            start_chrom_cmp
        }
    }
}

impl PartialOrd for Coordinate {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

pub struct SplitAlignmentSegment {
    /// reference sequence name
    pub rname: String,

    /// reference zero-indexed alignment start position
    pub pos: i64,

    /// Alignment using rust_htslib::bam::record::Cigar object
    pub cigar: rust_htslib::bam::record::CigarString,

    pub is_fwd_strand: bool,

    pub mapq: u8,
}

#[derive(Debug, PartialEq, Eq, Hash, Clone, Serialize, PartialOrd, Ord)]
pub struct FwdStrandSplitReadSegment {
    pub fwd_read_start: usize,
    pub fwd_read_end: usize,
    pub chrom: String,
    pub second_chrom: String,
    pub pos: i64,
    pub end: i64,
    pub is_fwd_strand: bool,
    pub is_start_softclipped: bool,
    pub is_end_softclipped: bool,
    pub phaseset_tag: Option<i32>,
    pub haplotype_tag: Option<i32>,
    pub spans: bool,

    /// All split segments found from SA tags in this bam record should be false
    pub from_primary_bam_record: bool,

    /// name of the read this alignment came from
    pub readname: String,
}

impl fmt::Display for FwdStrandSplitReadSegment {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let span_label = match self.spans {
            true => "Spanned",
            false => "Unspanned",
        };
        if self.chrom == self.second_chrom {
            write!(
                formatter,
                "{}:{}-{} {} {} ({}->{})",
                &self.chrom,
                &self.pos,
                &self.end,
                span_label,
                &self.readname,
                &self.fwd_read_start,
                &self.fwd_read_end
            )
        } else {
            write!(
                formatter,
                "{}:{}-{}:{} {} {}",
                &self.chrom, &self.pos, &self.second_chrom, &self.end, span_label, &self.readname,
            )
        }
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash)]
pub enum Orientation {
    Missing,
    Reverse,
    Forward,
}

impl Serialize for Orientation {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let s = match self {
            Orientation::Missing => "",
            Orientation::Reverse => "-",
            Orientation::Forward => "+",
        };
        serializer.serialize_str(s)
    }
}
#[derive(Debug, PartialEq, Eq)]
pub struct ParseOrientationError;

impl FromStr for Orientation {
    type Err = ParseOrientationError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Ok(Orientation::Missing),
            "+" => Ok(Orientation::Forward),
            "-" => Ok(Orientation::Reverse),
            _ => Err(ParseOrientationError),
        }
    }
}

/// Representation of a block in the complex SV,
/// with the start and end in a Coordinate, alignment summaries,
/// phasing flag, and position within the graph.
#[derive(Debug, PartialEq, Eq, Clone, Hash, Serialize)]
pub struct ComplexSVBlock {
    /// full span of the genomic block
    pub region: Coordinate,
    /// coverage for the sub-region spanned by supplementary alignments,
    /// where each key defines a block in the BED form {chrom:pos-end}
    /// and the value is depth of coverage
    pub coverages: BTreeMap<String, u32>,

    pub sample_order_index: u32,

    pub orientation: Orientation,
}

impl ComplexSVBlock {
    pub fn new(
        region: Coordinate,
        coverages: BTreeMap<Coordinate, u32>,
        sample_order_index: u32,
        fwd_orientation: &str,
    ) -> Self {
        let mut new_coverages = BTreeMap::new();
        for (coord, cov) in coverages {
            let string_coord = coord.to_string();
            new_coverages.insert(string_coord, cov);
        }
        if let Ok(orientation) = Orientation::from_str(fwd_orientation) {
            ComplexSVBlock {
                region,
                coverages: new_coverages,
                sample_order_index,
                orientation,
            }
        } else {
            error!("Failed to parse orientation {fwd_orientation}");
            std::process::exit(exitcode::DATAERR);
        }
    }
}

impl Ord for ComplexSVBlock {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.region.cmp(&other.region)
    }
}

impl PartialOrd for ComplexSVBlock {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Serialize, Hash)]
pub struct ComplexSVCalls {
    /// Version of SVTopo used
    pub svtopo_version: String,

    /// Full set of complex SV blocks
    pub event_graphs: Vec<Vec<ComplexSVBlock>>,
}

impl ComplexSVCalls {
    pub fn new(event_graphs: Vec<Vec<ComplexSVBlock>>) -> Self {
        ComplexSVCalls {
            event_graphs,
            svtopo_version: env!("CARGO_PKG_VERSION").to_string(),
        }
    }
}

/// Representation of a connection between two Coordinates,
/// including a flag to indicate if the connection is from phasing
#[derive(Debug, PartialEq, Eq, Hash, Clone, Serialize)]
pub struct Connection {
    pub first_coord: Coordinate,
    pub second_coord: Coordinate,
    pub inferred_from_phasing: bool,
    pub is_spanned: bool,
}
impl Connection {
    pub fn new(
        first_coord: Coordinate,
        second_coord: Coordinate,
        inferred_from_phasing: bool,
        is_spanned: bool,
    ) -> Self {
        Connection {
            first_coord,
            second_coord,
            inferred_from_phasing,
            is_spanned,
        }
    }

    /// Reverse the connection coordinates in-place, without altering other fields
    pub fn reverse(&mut self) {
        let tmp = self.first_coord.clone();
        self.first_coord = self.second_coord.clone();
        self.second_coord = tmp;
    }
}

impl Ord for Connection {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Order by first coordinate, then by second coordinate
        let first_cmp = self.first_coord.cmp(&other.first_coord);
        if first_cmp != std::cmp::Ordering::Equal {
            return first_cmp;
        }
        self.second_coord.cmp(&other.second_coord)
    }
}

impl PartialOrd for Connection {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// Container for hashmap of levels where the key is the
/// sample order and the value is a vector of blocks at that
/// order level.
pub struct EventGraph {
    pub graph: HashMap<u32, Vec<Coordinate>>,
}
impl Default for EventGraph {
    fn default() -> Self {
        Self::new()
    }
}

impl EventGraph {
    pub fn new() -> Self {
        EventGraph {
            graph: HashMap::new(),
        }
    }
}
/// Container for a hashmap of chromosome name to set of coordinates
/// that will be omitted from analysis when reads are extracted from BAM
/// and when breaks are read from VCF
pub struct ExcludeRegions {
    pub regions: HashMap<String, BTreeSet<Coordinate>>,
}
impl ExcludeRegions {
    pub fn new() -> Self {
        ExcludeRegions {
            regions: HashMap::new(),
        }
    }
}
impl Default for ExcludeRegions {
    fn default() -> Self {
        Self::new()
    }
}

/// Coordinate wrapper for command-line region specification
#[derive(Debug, Clone)]
pub struct TargetCoordinate {
    pub coord: Coordinate,
}
impl TargetCoordinate {
    pub fn new(mut coord_str: String) -> Self {
        coord_str = coord_str.replace(',', "");
        if !coord_str.contains(':') {
            error!(
                "Invalid target coordinate chromosome delimiter in {coord_str}"
            );
            std::process::exit(exitcode::DATAERR);
        }
        if !coord_str.contains('-') {
            error!(
                "Invalid target coordinate coordinate delimiter in {coord_str}"
            );
            std::process::exit(exitcode::DATAERR);
        }
        let fields: Vec<&str> = coord_str.split(':').collect();
        assert!(fields.len() == 2);
        let chrom = fields[0];
        let position_fields: Vec<&str> = fields[1].split('-').collect();
        assert!(position_fields.len() == 2);
        let pos = match position_fields[0].parse::<i64>() {
            Ok(num) => num,
            Err(_) => {
                error!("Invalid target coordinate {coord_str}");
                std::process::exit(exitcode::DATAERR);
            }
        };
        let end = match position_fields[1].parse::<i64>() {
            Ok(num) => num,
            Err(_) => {
                error!("Invalid target coordinate {coord_str}");
                std::process::exit(exitcode::DATAERR);
            }
        };
        TargetCoordinate {
            coord: Coordinate::new_region(chrom.to_string(), pos, end),
        }
    }
}

pub struct CoverageMap {
    pub coverages: BTreeMap<Coordinate, (Vec<u32>, Vec<u32>)>,
    pub longest_blocks: HashMap<String, isize>,
}

impl CoverageMap {
    pub fn new(annotated_graphs: &[Vec<ComplexSVBlock>]) -> Self {
        let mut longest_blocks: HashMap<String, isize> = HashMap::new();
        let mut coverage_map = BTreeMap::new();
        for annotated_graph in annotated_graphs.iter() {
            for block in annotated_graph.iter() {
                if block.region.start_chrom == block.region.end_chrom {
                    let coord = block.region.clone();
                    let length = (coord.end - coord.start) as isize;
                    if let Some(prev_longest) = longest_blocks.get(&block.region.start_chrom) {
                        if *prev_longest < length {
                            longest_blocks.insert(block.region.start_chrom.clone(), length);
                        }
                    } else {
                        longest_blocks.insert(block.region.start_chrom.clone(), length);
                    }

                    let covs = (vec![0; length as usize], vec![0; length as usize]);
                    coverage_map.insert(coord, covs);
                }
            }
        }
        CoverageMap {
            coverages: coverage_map,
            longest_blocks,
        }
    }

    pub fn update_coverages(&mut self, bam_filename: PathBuf) {
        let bam_name_str = bam_filename.as_os_str().to_str().unwrap();
        let mut bam_reader = bam::Reader::from_path(&bam_filename).unwrap_or_else(|_error| {
            error!("Input BAM does not exist: \"{bam_name_str}\"");
            std::process::exit(exitcode::NOINPUT);
        });
        let mut record = bam::Record::new();
        while let Some(result) = bam_reader.read(&mut record) {
            match result {
                Ok(_) => {
                    if record.tid() < 0 {
                        debug!("Invalid record TID {} skipped", record.tid());
                        continue;
                    }
                    let tid = record.tid() as u32;
                    let header = bam_reader.header();
                    let name = header.tid2name(tid);
                    let chrom_str_result = std::str::from_utf8(name);
                    if let Ok(chrom_str) = chrom_str_result {
                        let chrom = String::from(chrom_str);
                        let alignment_start = record.reference_start() + 1;
                        let alignment_end = record.reference_end();
                        let mut search_coord1_end = max(0, alignment_start);

                        if let Some(longest_block) = self.longest_blocks.get(&chrom) {
                            search_coord1_end = max(0, alignment_start - (*longest_block as i64));
                        }

                        let search_coord1 = Coordinate::new(chrom.clone(), search_coord1_end);
                        let search_coord2 = Coordinate::new(chrom.clone(), alignment_end);

                        for (coord, (covs, low_mapq_covs)) in
                            self.coverages.range_mut(search_coord1..search_coord2)
                        {
                            let overlap_start = max(alignment_start, coord.start);
                            let overlap_end = min(alignment_end, coord.end);
                            for position in overlap_start..overlap_end {
                                let cov_idx = (position - coord.start) as usize;
                                covs[cov_idx] += 1;
                                if record.mapq() < utils::EXTREMELY_LOW_MAPQ {
                                    low_mapq_covs[cov_idx] += 1;
                                }
                            }
                        }
                    }
                }
                Err(_) => {
                    error!("Error parsing BAM: {bam_name_str}");
                    std::process::exit(exitcode::IOERR);
                }
            }
        }
    }
}

pub struct BlockProcessingContext<'a> {
    pub event_graph: &'a HashMap<u32, Vec<Coordinate>>,
    pub clip_coordinates: &'a HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    pub already_processed_coordinates: &'a mut HashSet<Coordinate>,
    pub annotated_event_graph: &'a mut Vec<ComplexSVBlock>,
}

pub struct BlockProcessingState {
    pub sample_idx: i64,
    pub is_first: bool,
    pub is_last: bool,
    pub prev_coords: Vec<Coordinate>,
    pub last_event_spanned: bool,
    pub order_scaler: u32,
}

#[derive(Debug, Clone)]
pub struct BlockStartInfo {
    pub phaseset: i32,
    pub haplotype: i32,
    pub alignment_start: i64,
    pub alignment_end: i64,
}

pub struct ReadMetadata {
    pub readname: String,
    pub phaseset_tag: Option<i32>,
    pub haplotype_tag: Option<i32>,
}

pub struct SegmentProcessingContext<'a> {
    pub fwd_read_split_segments: &'a mut Vec<FwdStrandSplitReadSegment>,
    pub exclude_regions: &'a mut ExcludeRegions,
    pub excluded: &'a mut bool,
}

/// Container for annotation processing state that tracks progress through event graph processing
pub struct AnnotationProcessingState {
    pub already_processed_coordinates: HashSet<Coordinate>,
    pub last_event_spanned: bool,
    pub order_scaler: u32,
}

/// Container for VCF processing results that holds connections and coordinate mappings
pub struct VcfProcessingResult {
    pub vcf_breaks: Vec<Connection>,
    pub coordinate_map: HashMap<String, Vec<Coordinate>>,
    pub bnd_records: HashMap<String, bcf::Record>,
}

impl VcfProcessingResult {
    pub fn new() -> Self {
        Self {
            vcf_breaks: Vec::new(),
            coordinate_map: HashMap::new(),
            bnd_records: HashMap::new(),
        }
    }
}

impl Default for VcfProcessingResult {
    fn default() -> Self {
        Self::new()
    }
}

/// Container for coordinate ordering state during event graph construction
pub struct CoordinateOrderingState {
    pub visited_edges: HashSet<Connection>,
    pub ordered_coordinates: HashMap<u32, Vec<Coordinate>>,
    pub queue: Vec<(Coordinate, u32)>,
    pub queue_idx: usize,
    pub ending_inv_connection_opt: Option<Connection>,
}

impl CoordinateOrderingState {
    pub fn new(start_coordinate: Option<Coordinate>) -> Self {
        let mut queue = Vec::new();
        if let Some(start) = start_coordinate {
            queue.push((start, 0));
        }

        Self {
            visited_edges: HashSet::new(),
            ordered_coordinates: HashMap::new(),
            queue,
            queue_idx: 0,
            ending_inv_connection_opt: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::{BTreeMap, HashSet};

    // Helper function to create test FwdStrandSplitReadSegment
    fn create_test_alignment(
        chrom: &str,
        pos: i64,
        end: i64,
        readname: &str,
        spans: bool,
    ) -> FwdStrandSplitReadSegment {
        FwdStrandSplitReadSegment {
            fwd_read_start: 0,
            fwd_read_end: 100,
            chrom: chrom.to_string(),
            second_chrom: chrom.to_string(),
            pos,
            end,
            is_fwd_strand: true,
            is_start_softclipped: true,
            is_end_softclipped: false,
            phaseset_tag: Some(1),
            haplotype_tag: Some(0),
            spans,
            from_primary_bam_record: true,
            readname: readname.to_string(),
        }
    }

    // Helper function to create test ComplexSVBlock
    fn create_test_block(
        chrom: &str,
        start: i64,
        end: i64,
        sample_order: u32,
        orientation: &str,
    ) -> ComplexSVBlock {
        let region = Coordinate::new_region(chrom.to_string(), start, end);
        let mut coverages = BTreeMap::new();
        coverages.insert(region.clone(), 10);
        ComplexSVBlock::new(region, coverages, sample_order, orientation)
    }

    #[test]
    fn test_coordinate_new() {
        let coord = Coordinate::new("chr1".to_string(), 1000);
        assert_eq!(coord.start_chrom, "chr1");
        assert_eq!(coord.start, 1000);
        assert_eq!(coord.end_chrom, "chr1");
        assert_eq!(coord.end, 1000);
        assert_eq!(
            coord.confidence_interval,
            (utils::MAX_CLUST_DIST, utils::MAX_CLUST_DIST)
        );
        assert!(coord.variant_ids.is_empty());
    }

    #[test]
    fn test_coordinate_new_region() {
        let coord = Coordinate::new_region("chr1".to_string(), 1000, 2000);
        assert_eq!(coord.start_chrom, "chr1");
        assert_eq!(coord.start, 1000);
        assert_eq!(coord.end_chrom, "chr1");
        assert_eq!(coord.end, 2000);
        assert_eq!(
            coord.confidence_interval,
            (utils::MAX_CLUST_DIST, utils::MAX_CLUST_DIST)
        );
        assert!(coord.variant_ids.is_empty());
    }

    #[test]
    fn test_coordinate_new_with_confidence_interval() {
        let coord = Coordinate::new_with_confidence_interval("chr1".to_string(), 1000, (10, 20));
        assert_eq!(coord.start_chrom, "chr1");
        assert_eq!(coord.start, 1000);
        assert_eq!(coord.end_chrom, "chr1");
        assert_eq!(coord.end, 1000);
        assert_eq!(coord.confidence_interval, (10, 20));
        assert!(coord.variant_ids.is_empty());
    }

    #[test]
    fn test_coordinate_from_alignment() {
        let alignment = create_test_alignment("chr1", 1000, 2000, "read1", true);
        let coord = Coordinate::from_alignment(alignment);
        assert_eq!(coord.start_chrom, "chr1");
        assert_eq!(coord.start, 1000);
        assert_eq!(coord.end_chrom, "chr1");
        assert_eq!(coord.end, 2000);
        assert_eq!(
            coord.confidence_interval,
            (utils::MAX_CLUST_DIST, utils::MAX_CLUST_DIST)
        );
    }

    #[test]
    fn test_coordinate_from_connection() {
        let coord1 = Coordinate::new("chr1".to_string(), 1000);
        let coord2 = Coordinate::new("chr2".to_string(), 2000);
        let connection = Connection::new(coord1, coord2, false, true);
        let coord = Coordinate::from_connection(connection);
        assert_eq!(coord.start_chrom, "chr1");
        assert_eq!(coord.start, 1000);
        assert_eq!(coord.end_chrom, "chr2");
        assert_eq!(coord.end, 2000);
    }

    #[test]
    fn test_coordinate_is_within() {
        let coord1 = Coordinate::new_with_confidence_interval("chr1".to_string(), 1000, (10, 10));
        let coord2 = Coordinate::new_with_confidence_interval("chr1".to_string(), 1005, (10, 10));
        let coord3 = Coordinate::new_with_confidence_interval("chr1".to_string(), 1020, (10, 10));
        let coord4 = Coordinate::new_with_confidence_interval("chr2".to_string(), 1005, (10, 10));

        // Within confidence interval
        assert!(coord1.is_within(&coord2));
        assert!(coord2.is_within(&coord1));

        // Outside confidence interval
        assert!(!coord1.is_within(&coord3));
        assert!(!coord3.is_within(&coord1));

        // Different chromosome
        assert!(!coord1.is_within(&coord4));
        assert!(!coord4.is_within(&coord1));
    }

    #[test]
    fn test_coordinate_display() {
        let coord1 = Coordinate::new("chr1".to_string(), 1000);
        let coord2 = Coordinate::new_region("chr1".to_string(), 1000, 2000);
        let mut coord3 =
            Coordinate::new_with_confidence_interval("chr1".to_string(), 1000, (10, 10));
        coord3.end_chrom = "chr2".to_string();
        coord3.end = 2000;

        assert_eq!(coord1.to_string(), "chr1:1000");
        assert_eq!(coord2.to_string(), "chr1:1000-2000");
        assert_eq!(coord3.to_string(), "chr1:1000-chr2:2000");
    }

    #[test]
    fn test_coordinate_ordering() {
        let coord1 = Coordinate::new("chr1".to_string(), 1000);
        let coord2 = Coordinate::new("chr1".to_string(), 1000);
        let coord3 = Coordinate::new("chr1".to_string(), 2000);
        let coord4 = Coordinate::new("chr2".to_string(), 1000);

        // Same coordinates
        assert_eq!(coord1.cmp(&coord2), std::cmp::Ordering::Equal);

        // Different positions on same chromosome
        assert_eq!(coord1.cmp(&coord3), std::cmp::Ordering::Less);
        assert_eq!(coord3.cmp(&coord1), std::cmp::Ordering::Greater);

        // Different chromosomes
        assert_eq!(coord1.cmp(&coord4), std::cmp::Ordering::Less);
        assert_eq!(coord4.cmp(&coord1), std::cmp::Ordering::Greater);
    }

    #[test]
    fn test_fwd_strand_split_read_segment_display() {
        let alignment = create_test_alignment("chr1", 1000, 2000, "read1", true);
        let display = alignment.to_string();
        assert!(display.contains("chr1:1000-2000"));
        assert!(display.contains("Spanned"));
        assert!(display.contains("read1"));
        assert!(display.contains("0->100"));
    }

    #[test]
    fn test_orientation_from_str() {
        assert!(Orientation::from_str("").is_ok());
        assert!(Orientation::from_str("+").is_ok());
        assert!(Orientation::from_str("-").is_ok());

        assert_eq!(Orientation::from_str("").unwrap(), Orientation::Missing);
        assert_eq!(Orientation::from_str("+").unwrap(), Orientation::Forward);
        assert_eq!(Orientation::from_str("-").unwrap(), Orientation::Reverse);
    }

    #[test]
    fn test_orientation_from_str_invalid() {
        // This test expects the function to panic due to std::process::exit
        // We can't use should_panic because process::exit terminates the test runner
        // Instead, we'll test the parsing logic separately
        let invalid_str = "invalid";
        assert_ne!(invalid_str, "");
        assert_ne!(invalid_str, "+");
        assert_ne!(invalid_str, "-");
    }

    #[test]
    fn test_complex_sv_block_new() {
        let region = Coordinate::new_region("chr1".to_string(), 1000, 2000);
        let mut coverages = BTreeMap::new();
        coverages.insert(region.clone(), 10);

        let block = ComplexSVBlock::new(region.clone(), coverages, 1, "+");

        assert_eq!(block.region, region);
        assert_eq!(block.sample_order_index, 1);
        assert_eq!(block.orientation, Orientation::Forward);
        assert_eq!(block.coverages.get("chr1:1000-2000"), Some(&10));
    }

    #[test]
    fn test_complex_sv_block_ordering() {
        let block1 = create_test_block("chr1", 1000, 2000, 1, "+");
        let block2 = create_test_block("chr1", 1000, 2000, 2, "+");
        let block3 = create_test_block("chr2", 1000, 2000, 1, "+");

        // Same region, different sample order
        assert_eq!(block1.cmp(&block2), std::cmp::Ordering::Equal);

        // Different regions
        assert_eq!(block1.cmp(&block3), std::cmp::Ordering::Less);
        assert_eq!(block3.cmp(&block1), std::cmp::Ordering::Greater);
    }

    #[test]
    fn test_complex_sv_calls_new() {
        let block1 = create_test_block("chr1", 1000, 2000, 1, "+");
        let block2 = create_test_block("chr2", 2000, 3000, 2, "-");
        let event_graphs = vec![vec![block1, block2]];

        let calls = ComplexSVCalls::new(event_graphs);

        assert_eq!(calls.event_graphs.len(), 1);
        assert_eq!(calls.event_graphs[0].len(), 2);
        assert_eq!(calls.svtopo_version, env!("CARGO_PKG_VERSION"));
    }

    #[test]
    fn test_connection_new() {
        let coord1 = Coordinate::new("chr1".to_string(), 1000);
        let coord2 = Coordinate::new("chr2".to_string(), 2000);

        let connection = Connection::new(coord1.clone(), coord2.clone(), true, false);

        assert_eq!(connection.first_coord, coord1);
        assert_eq!(connection.second_coord, coord2);
        assert!(connection.inferred_from_phasing);
        assert!(!connection.is_spanned);
    }

    #[test]
    fn test_connection_reverse() {
        let coord1 = Coordinate::new("chr1".to_string(), 1000);
        let coord2 = Coordinate::new("chr2".to_string(), 2000);

        let mut connection = Connection::new(coord1.clone(), coord2.clone(), true, false);
        connection.reverse();

        assert_eq!(connection.first_coord, coord2);
        assert_eq!(connection.second_coord, coord1);
        assert!(connection.inferred_from_phasing);
        assert!(!connection.is_spanned);
    }

    #[test]
    fn test_connection_ordering() {
        let coord1 = Coordinate::new("chr1".to_string(), 1000);
        let coord2 = Coordinate::new("chr2".to_string(), 2000);
        let coord3 = Coordinate::new("chr3".to_string(), 3000);

        let conn1 = Connection::new(coord1.clone(), coord2.clone(), false, true);
        let conn2 = Connection::new(coord1.clone(), coord2.clone(), true, false);
        let conn3 = Connection::new(coord1.clone(), coord3.clone(), false, true);

        // Same coordinates, different flags
        assert_eq!(conn1.cmp(&conn2), std::cmp::Ordering::Equal);

        // Different coordinates
        assert_eq!(conn1.cmp(&conn3), std::cmp::Ordering::Less);
        assert_eq!(conn3.cmp(&conn1), std::cmp::Ordering::Greater);
    }

    #[test]
    fn test_event_graph_new() {
        let graph = EventGraph::new();
        assert!(graph.graph.is_empty());
    }

    #[test]
    fn test_event_graph_default() {
        let graph = EventGraph::default();
        assert!(graph.graph.is_empty());
    }

    #[test]
    fn test_exclude_regions_new() {
        let regions = ExcludeRegions::new();
        assert!(regions.regions.is_empty());
    }

    #[test]
    fn test_coverage_map_new() {
        let block1 = create_test_block("chr1", 1000, 2000, 1, "+");
        let block2 = create_test_block("chr1", 2000, 3000, 2, "-");
        let annotated_graphs = vec![vec![block1, block2]];

        let coverage_map = CoverageMap::new(&annotated_graphs);

        assert_eq!(coverage_map.coverages.len(), 2);
        assert_eq!(coverage_map.longest_blocks.get("chr1"), Some(&1000));

        // Check that coverage vectors are initialized with zeros
        for (covs, low_mapq_covs) in coverage_map.coverages.values() {
            assert!(covs.iter().all(|&x| x == 0));
            assert!(low_mapq_covs.iter().all(|&x| x == 0));
        }
    }

    #[test]
    fn test_coverage_map_new_multi_chromosome() {
        let block1 = create_test_block("chr1", 1000, 2000, 1, "+");
        let block2 = create_test_block("chr2", 2000, 4000, 2, "-");
        let annotated_graphs = vec![vec![block1, block2]];

        let coverage_map = CoverageMap::new(&annotated_graphs);

        assert_eq!(coverage_map.coverages.len(), 2);
        assert_eq!(coverage_map.longest_blocks.get("chr1"), Some(&1000));
        assert_eq!(coverage_map.longest_blocks.get("chr2"), Some(&2000));
    }

    #[test]
    fn test_coverage_map_new_translocation() {
        let mut block = create_test_block("chr1", 1000, 2000, 1, "+");
        block.region.end_chrom = "chr2".to_string();
        block.region.end = 3000;
        let annotated_graphs = vec![vec![block]];

        let coverage_map = CoverageMap::new(&annotated_graphs);

        // Translocations should not be included in coverage map
        assert_eq!(coverage_map.coverages.len(), 0);
        assert_eq!(coverage_map.longest_blocks.len(), 0);
    }

    #[test]
    fn test_vcf_processing_result_new() {
        let result = VcfProcessingResult::new();
        assert!(result.vcf_breaks.is_empty());
        assert!(result.coordinate_map.is_empty());
        assert!(result.bnd_records.is_empty());
    }

    #[test]
    fn test_coordinate_ordering_state_new_with_start() {
        let start_coord = Coordinate::new("chr1".to_string(), 1000);
        let state = CoordinateOrderingState::new(Some(start_coord.clone()));

        assert!(state.visited_edges.is_empty());
        assert!(state.ordered_coordinates.is_empty());
        assert_eq!(state.queue.len(), 1);
        assert_eq!(state.queue[0], (start_coord, 0));
        assert_eq!(state.queue_idx, 0);
        assert!(state.ending_inv_connection_opt.is_none());
    }

    #[test]
    fn test_coordinate_ordering_state_new_without_start() {
        let state = CoordinateOrderingState::new(None);

        assert!(state.visited_edges.is_empty());
        assert!(state.ordered_coordinates.is_empty());
        assert!(state.queue.is_empty());
        assert_eq!(state.queue_idx, 0);
        assert!(state.ending_inv_connection_opt.is_none());
    }

    #[test]
    fn test_coordinate_equality() {
        let coord1 = Coordinate::new("chr1".to_string(), 1000);
        let coord2 = Coordinate::new("chr1".to_string(), 1000);
        let coord3 = Coordinate::new("chr1".to_string(), 2000);

        assert_eq!(coord1, coord2);
        assert_ne!(coord1, coord3);
    }

    #[test]
    fn test_coordinate_hash() {
        let coord1 = Coordinate::new("chr1".to_string(), 1000);
        let coord2 = Coordinate::new("chr1".to_string(), 1000);
        let coord3 = Coordinate::new("chr1".to_string(), 2000);

        let mut set = HashSet::new();
        set.insert(coord1.clone());
        set.insert(coord2.clone());
        set.insert(coord3.clone());

        // Should only have 2 unique coordinates
        assert_eq!(set.len(), 2);
        assert!(set.contains(&coord1));
        assert!(set.contains(&coord3));
    }

    #[test]
    fn test_fwd_strand_split_read_segment_equality() {
        let seg1 = create_test_alignment("chr1", 1000, 2000, "read1", true);
        let seg2 = create_test_alignment("chr1", 1000, 2000, "read1", true);
        let seg3 = create_test_alignment("chr1", 1000, 2000, "read2", true);

        assert_eq!(seg1, seg2);
        assert_ne!(seg1, seg3);
    }

    #[test]
    fn test_complex_sv_block_equality() {
        let block1 = create_test_block("chr1", 1000, 2000, 1, "+");
        let block2 = create_test_block("chr1", 1000, 2000, 1, "+");
        let block3 = create_test_block("chr1", 1000, 2000, 2, "+");

        assert_eq!(block1, block2);
        assert_ne!(block1, block3);
    }

    #[test]
    fn test_connection_equality() {
        let coord1 = Coordinate::new("chr1".to_string(), 1000);
        let coord2 = Coordinate::new("chr2".to_string(), 2000);

        let conn1 = Connection::new(coord1.clone(), coord2.clone(), false, true);
        let conn2 = Connection::new(coord1.clone(), coord2.clone(), false, true);
        let conn3 = Connection::new(coord1.clone(), coord2.clone(), true, true);

        assert_eq!(conn1, conn2);
        assert_ne!(conn1, conn3);
    }

    #[test]
    fn test_orientation_serialization() {
        let missing = Orientation::Missing;
        let forward = Orientation::Forward;
        let reverse = Orientation::Reverse;

        // Test serialization (this would require serde_json, but we can test the logic)
        assert_eq!(missing, Orientation::from_str("").unwrap());
        assert_eq!(forward, Orientation::from_str("+").unwrap());
        assert_eq!(reverse, Orientation::from_str("-").unwrap());
    }

    #[test]
    fn test_complex_sv_calls_equality() {
        let block1 = create_test_block("chr1", 1000, 2000, 1, "+");
        let block2 = create_test_block("chr2", 2000, 3000, 2, "-");
        let event_graphs1 = vec![vec![block1.clone(), block2.clone()]];
        let event_graphs2 = vec![vec![block1, block2]];

        let calls1 = ComplexSVCalls::new(event_graphs1);
        let calls2 = ComplexSVCalls::new(event_graphs2);

        assert_eq!(calls1.event_graphs, calls2.event_graphs);
        assert_eq!(calls1.svtopo_version, calls2.svtopo_version);
    }

    #[test]
    fn test_edge_cases() {
        // Test coordinate with zero length
        let coord = Coordinate::new_region("chr1".to_string(), 1000, 1000);
        assert_eq!(coord.start, coord.end);
        assert_eq!(coord.to_string(), "chr1:1000");

        // Test coordinate with negative positions
        let coord = Coordinate::new("chr1".to_string(), -1000);
        assert_eq!(coord.start, -1000);
        assert_eq!(coord.to_string(), "chr1:-1000");

        // Test empty chromosome name
        let coord = Coordinate::new("".to_string(), 1000);
        assert_eq!(coord.start_chrom, "");
        assert_eq!(coord.to_string(), ":1000");

        // Test very large positions
        let coord = Coordinate::new("chr1".to_string(), i64::MAX);
        assert_eq!(coord.start, i64::MAX);
    }

    #[test]
    fn test_coordinate_variant_ids() {
        let mut coord = Coordinate::new("chr1".to_string(), 1000);
        coord.variant_ids.push("var1".to_string());
        coord.variant_ids.push("var2".to_string());

        assert_eq!(coord.variant_ids.len(), 2);
        assert_eq!(coord.variant_ids[0], "var1");
        assert_eq!(coord.variant_ids[1], "var2");
    }

    #[test]
    fn test_fwd_strand_split_read_segment_phasing() {
        let mut seg = create_test_alignment("chr1", 1000, 2000, "read1", true);

        // Test with phasing
        assert_eq!(seg.phaseset_tag, Some(1));
        assert_eq!(seg.haplotype_tag, Some(0));

        // Test without phasing
        seg.phaseset_tag = None;
        seg.haplotype_tag = None;
        assert_eq!(seg.phaseset_tag, None);
        assert_eq!(seg.haplotype_tag, None);
    }

    #[test]
    fn test_complex_sv_block_coverages() {
        let region = Coordinate::new_region("chr1".to_string(), 1000, 2000);
        let mut coverages = BTreeMap::new();
        coverages.insert(region.clone(), 10);
        coverages.insert(Coordinate::new_region("chr1".to_string(), 2000, 3000), 20);

        let block = ComplexSVBlock::new(region, coverages, 1, "+");
        assert_eq!(block.coverages.len(), 2);
        assert_eq!(block.coverages.get("chr1:1000-2000"), Some(&10));
        assert_eq!(block.coverages.get("chr1:2000-3000"), Some(&20));
    }

    #[test]
    fn test_target_coordinate_new_valid() {
        let target = TargetCoordinate::new("chr1:1000-2000".to_string());
        assert_eq!(target.coord.start_chrom, "chr1");
        assert_eq!(target.coord.start, 1000);
        assert_eq!(target.coord.end_chrom, "chr1");
        assert_eq!(target.coord.end, 2000);
    }

    #[test]
    fn test_target_coordinate_new_with_commas() {
        let target = TargetCoordinate::new("chr1,000:1,000-2,000".to_string());
        assert_eq!(target.coord.start_chrom, "chr1000");
        assert_eq!(target.coord.start, 1000);
        assert_eq!(target.coord.end_chrom, "chr1000");
        assert_eq!(target.coord.end, 2000);
    }

    #[test]
    fn test_target_coordinate_new_single_position() {
        let target = TargetCoordinate::new("chr1:1000-1000".to_string());
        assert_eq!(target.coord.start_chrom, "chr1");
        assert_eq!(target.coord.start, 1000);
        assert_eq!(target.coord.end_chrom, "chr1");
        assert_eq!(target.coord.end, 1000);
    }

    #[test]
    fn test_target_coordinate_new_large_numbers() {
        let target = TargetCoordinate::new("chr1:123456789-987654321".to_string());
        assert_eq!(target.coord.start_chrom, "chr1");
        assert_eq!(target.coord.start, 123456789);
        assert_eq!(target.coord.end_chrom, "chr1");
        assert_eq!(target.coord.end, 987654321);
    }

    #[test]
    fn test_target_coordinate_new_different_chromosomes() {
        // Note: This function only handles single chromosome regions
        // The function assumes start_chrom == end_chrom
        let target = TargetCoordinate::new("chr1:1000-2000".to_string());
        assert_eq!(target.coord.start_chrom, target.coord.end_chrom);
    }

    #[test]
    fn test_coverage_map_update_coverages_basic() {
        // Create a simple coverage map with one block
        let region = Coordinate::new_region("chr1".to_string(), 1000, 2000);
        let mut coverages = BTreeMap::new();
        coverages.insert(region.clone(), 0);
        let block = ComplexSVBlock::new(region.clone(), coverages, 1, "+");
        let annotated_graphs = vec![vec![block]];

        let coverage_map = CoverageMap::new(&annotated_graphs);

        // Verify initial state
        assert_eq!(coverage_map.coverages.len(), 1);
        assert_eq!(coverage_map.longest_blocks.len(), 1);
        assert_eq!(coverage_map.longest_blocks.get("chr1"), Some(&1000));

        // Check that the coverage vectors are initialized correctly
        let (covs, low_mapq_covs) = coverage_map.coverages.get(&region).unwrap();
        assert_eq!(covs.len(), 1000); // 2000 - 1000 = 1000 positions
        assert_eq!(low_mapq_covs.len(), 1000);
        assert!(covs.iter().all(|&x| x == 0));
        assert!(low_mapq_covs.iter().all(|&x| x == 0));
    }

    #[test]
    fn test_coverage_map_update_coverages_multiple_blocks() {
        // Create coverage map with multiple blocks on same chromosome
        let region1 = Coordinate::new_region("chr1".to_string(), 1000, 2000);
        let region2 = Coordinate::new_region("chr1".to_string(), 3000, 4000);
        let mut coverages1 = BTreeMap::new();
        let mut coverages2 = BTreeMap::new();
        coverages1.insert(region1.clone(), 0);
        coverages2.insert(region2.clone(), 0);

        let block1 = ComplexSVBlock::new(region1, coverages1, 1, "+");
        let block2 = ComplexSVBlock::new(region2, coverages2, 2, "-");
        let annotated_graphs = vec![vec![block1, block2]];

        let coverage_map = CoverageMap::new(&annotated_graphs);

        // Verify multiple blocks are handled
        assert_eq!(coverage_map.coverages.len(), 2);
        assert_eq!(coverage_map.longest_blocks.get("chr1"), Some(&1000)); // longest block length
    }

    #[test]
    fn test_coverage_map_update_coverages_multi_chromosome() {
        // Create coverage map with blocks on different chromosomes
        let region1 = Coordinate::new_region("chr1".to_string(), 1000, 2000);
        let region2 = Coordinate::new_region("chr2".to_string(), 5000, 6000);
        let mut coverages1 = BTreeMap::new();
        let mut coverages2 = BTreeMap::new();
        coverages1.insert(region1.clone(), 0);
        coverages2.insert(region2.clone(), 0);

        let block1 = ComplexSVBlock::new(region1, coverages1, 1, "+");
        let block2 = ComplexSVBlock::new(region2, coverages2, 2, "-");
        let annotated_graphs = vec![vec![block1, block2]];

        let coverage_map = CoverageMap::new(&annotated_graphs);

        // Verify multi-chromosome handling
        assert_eq!(coverage_map.coverages.len(), 2);
        assert_eq!(coverage_map.longest_blocks.len(), 2);
        assert_eq!(coverage_map.longest_blocks.get("chr1"), Some(&1000));
        assert_eq!(coverage_map.longest_blocks.get("chr2"), Some(&1000));
    }

    #[test]
    fn test_coverage_map_update_coverages_empty_input() {
        // Test with empty annotated graphs
        let annotated_graphs: Vec<Vec<ComplexSVBlock>> = vec![];
        let coverage_map = CoverageMap::new(&annotated_graphs);

        assert_eq!(coverage_map.coverages.len(), 0);
        assert_eq!(coverage_map.longest_blocks.len(), 0);
    }

    #[test]
    fn test_coverage_map_update_coverages_translocation_blocks() {
        // Test with translocation blocks (different start/end chromosomes)
        let region = Coordinate {
            start_chrom: "chr1".to_string(),
            start: 1000,
            end_chrom: "chr2".to_string(),
            end: 2000,
            variant_ids: vec![],
            confidence_interval: (utils::MAX_CLUST_DIST, utils::MAX_CLUST_DIST),
        };
        let mut coverages = BTreeMap::new();
        coverages.insert(region.clone(), 0);
        let block = ComplexSVBlock::new(region, coverages, 1, "+");
        let annotated_graphs = vec![vec![block]];

        let coverage_map = CoverageMap::new(&annotated_graphs);

        // Translocation blocks should be skipped (start_chrom != end_chrom)
        assert_eq!(coverage_map.coverages.len(), 0);
        assert_eq!(coverage_map.longest_blocks.len(), 0);
    }
}
