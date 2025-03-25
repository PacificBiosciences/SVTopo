use core::fmt;
use log::{debug, error};
use rust_htslib::bam::{self, ext::BamRecordExtensions, Read};
use serde::Serialize;
use std::{
    cmp::{max, min},
    collections::{BTreeMap, BTreeSet, HashMap},
    path::PathBuf,
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

#[derive(Debug, PartialEq, Eq, Clone, Hash)]
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

impl Orientation {
    pub fn from_str(s: &str) -> Self {
        match s {
            "" => Orientation::Missing,
            "+" => Orientation::Forward,
            "-" => Orientation::Reverse,
            _ => {
                error!("Invalid value {} for Orientation", s);
                std::process::exit(exitcode::DATAERR);
            }
        }
    }
}

/// Representation of a block in the complex SV,
/// with the start and end in a Coordinate, alignment summaries,
/// phasing flag, and position within the graph.
#[derive(Debug, PartialEq, Eq, Clone, Serialize, Hash)]
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

        ComplexSVBlock {
            region,
            coverages: new_coverages,
            sample_order_index,
            orientation: Orientation::from_str(fwd_orientation),
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
#[derive(Debug, PartialEq, Eq, Hash, Clone, Serialize, PartialOrd, Ord)]
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
                "Invalid target coordinate chromosome delimiter in {}",
                coord_str
            );
            std::process::exit(exitcode::DATAERR);
        }
        if !coord_str.contains('-') {
            error!(
                "Invalid target coordinate coordinate delimiter in {}",
                coord_str
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
                error!("Invalid target coordinate {}", coord_str);
                std::process::exit(exitcode::DATAERR);
            }
        };
        let end = match position_fields[1].parse::<i64>() {
            Ok(num) => num,
            Err(_) => {
                error!("Invalid target coordinate {}", coord_str);
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
            error!("Input BAM does not exist: \"{}\"", bam_name_str);
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
                    error!("Error parsing BAM: {}", bam_name_str);
                    std::process::exit(exitcode::IOERR);
                }
            }
        }
    }
}
