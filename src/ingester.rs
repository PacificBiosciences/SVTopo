use flate2::bufread::GzDecoder;
use log::{debug, error};
use rust_htslib::bam;

use rust_htslib::bam::Read;
use rust_htslib::bcf;
use std::cmp::max;
use std::collections::HashMap;
use std::collections::HashSet;
use std::io::BufReader;
use std::path::Path;
use std::path::PathBuf;
use std::time::SystemTime;

use crate::bam_sa_parser::get_fwd_read_split_segments;
use crate::containers::TargetCoordinate;
use crate::containers::{
    Connection, Coordinate, ExcludeRegions, FwdStrandSplitReadSegment, VcfProcessingResult,
};
use crate::utils;
use crate::utils::is_gzipped;
use crate::utils::is_local_file;

/// Extract sample name from bam header
///
/// This uses the sample name from the first read group found in the header, and does not
/// check for additional read groups. Termination triggered if no sample ID found.
pub fn get_sample_from_bam(bam_filename: PathBuf) -> String {
    let bam_name_str = bam_filename.as_os_str().to_str().unwrap();
    let bam_reader: bam::Reader = bam::Reader::from_path(&bam_filename).unwrap_or_else(|_error| {
        error!("Input BAM does not exist: \"{bam_name_str}\"");
        std::process::exit(exitcode::NOINPUT);
    });

    for line in std::str::from_utf8(bam_reader.header().as_bytes())
        .unwrap()
        .split('\n')
    {
        for (i, word) in line.split('\t').enumerate() {
            if i == 0 {
                if word != "@RG" {
                    break;
                }
            } else if let Some(sample_name) = word.strip_prefix("SM:") {
                return sample_name.to_string();
            }
        }
    }

    error!("Input BAM {bam_name_str} missing sample name (SM tag)");
    std::process::exit(exitcode::DATAERR);
}

pub fn get_split_alignments_from_region(
    bam_filename: PathBuf,
    target_coordinate: TargetCoordinate,
) -> HashMap<String, Vec<FwdStrandSplitReadSegment>> {
    let start_time = SystemTime::now();
    let bam_name_str = bam_filename.as_os_str().to_str().unwrap();
    let mut bam_reader: bam::IndexedReader = bam::IndexedReader::from_path(&bam_filename)
        .unwrap_or_else(|_error| {
            error!("Input BAM does not exist: \"{bam_name_str}\"");
            std::process::exit(exitcode::NOINPUT);
        });

    let mut alignment_map: HashMap<String, Vec<FwdStrandSplitReadSegment>> = HashMap::new();
    // while more verbose, reading into an existing record rather than allocating new ones should
    // improve efficiency, since this method will only allocate memory for one record
    let fetch_region = (
        target_coordinate.coord.start_chrom.as_str(),
        target_coordinate.coord.start,
        target_coordinate.coord.end,
    );
    let _ = bam_reader.fetch(fetch_region);
    for record_result in bam_reader.records() {
        match record_result {
            Ok(record) => {
                if record.mapq() < utils::MIN_MAPQ {
                    continue;
                }
                // process only reads with significant softclipping
                if record.cigar().leading_softclips() < utils::MIN_CLIPPING
                    && record.cigar().trailing_softclips() < utils::MIN_CLIPPING
                {
                    continue;
                }
                let chrom = target_coordinate.coord.start_chrom.clone();
                let readname = String::from_utf8(record.qname().to_vec()).unwrap();
                let mut alignments =
                    get_fwd_read_split_segments(&record, chrom, &mut ExcludeRegions::new());
                if alignments.len() > 1 {
                    flip_reverse_alignments(&mut alignments);
                    alignment_map.insert(readname, alignments);
                }
            }
            Err(_) => {
                error!("Error parsing BAM: {bam_name_str}");
                std::process::exit(exitcode::IOERR);
            }
        }
    }
    debug!(
        "Read extraction: {}s",
        start_time.elapsed().unwrap().as_secs()
    );
    debug!("{} clipped reads", alignment_map.len());
    let alignment_count: usize = alignment_map.values().map(|s| s.len()).sum();
    debug!("{alignment_count} total chimeric alignments");
    alignment_map
}

/// Gets clipped reads where each alignment has at least 100 bases of softclipped
/// and there are at least two alignments.
///
/// This includes supplementary alignments and the single primary associated with each group of them
pub fn get_split_alignments(
    bam_filename: PathBuf,
    exclude_regions: &mut ExcludeRegions,
) -> HashMap<String, Vec<FwdStrandSplitReadSegment>> {
    let start_time = SystemTime::now();
    let bam_name_str = bam_filename.as_os_str().to_str().unwrap();
    let mut bam_reader: bam::Reader =
        bam::Reader::from_path(&bam_filename).unwrap_or_else(|_error| {
            error!("Input BAM does not exist: \"{bam_name_str}\"");
            std::process::exit(exitcode::NOINPUT);
        });

    let mut alignment_map: HashMap<String, Vec<FwdStrandSplitReadSegment>> = HashMap::new();
    // while more verbose, reading into an existing record rather than allocating new ones should
    // improve efficiency, since this method will only allocate memory for one record
    let mut record: bam::Record = bam::Record::new();
    while let Some(result) = bam::Read::read(&mut bam_reader, &mut record) {
        match result {
            Ok(_) => {
                if record.mapq() < utils::MIN_MAPQ {
                    continue;
                }
                // process only reads with significant softclipping
                if record.cigar().leading_softclips() < utils::MIN_CLIPPING
                    && record.cigar().trailing_softclips() < utils::MIN_CLIPPING
                {
                    continue;
                }
                let chrom = get_chrom_from_bam(&bam_reader, &record);
                let readname = String::from_utf8(record.qname().to_vec()).unwrap();
                let mut alignments = get_fwd_read_split_segments(&record, chrom, exclude_regions);
                if alignments.len() > 1 {
                    flip_reverse_alignments(&mut alignments);
                    alignment_map.insert(readname, alignments);
                }
            }
            Err(_) => {
                error!("Error parsing BAM: {bam_name_str}");
                std::process::exit(exitcode::IOERR);
            }
        }
    }
    debug!(
        "Read extraction: {}s",
        start_time.elapsed().unwrap().as_secs()
    );
    debug!("{} clipped reads", alignment_map.len());
    let alignment_count: usize = alignment_map.values().map(|s| s.len()).sum();
    debug!("{alignment_count} total chimeric alignments");
    alignment_map
}

/// get the chromosome name from the header, adding the chr
/// prefix if the genome is a known human reference.
fn get_chrom_from_bam(bam_reader: &bam::Reader, record: &bam::Record) -> String {
    let tid = record.tid();
    assert!(tid >= 0);
    let chrom_bytes = bam::Read::header(bam_reader).tid2name(tid as u32);
    let chrom = String::from_utf8(chrom_bytes.to_vec()).unwrap();

    // if the genome reference is unrecognized, we can't assume it should or shouldn't have a chr prefix
    // and if the reference is known but that prefix is already present, we need make no change
    if chrom.get(0..3) == Some("chr") {
        chrom
    } else {
        format!("chr{chrom}")
    }
}

/// Stores any connections that are found, whether from mated BNDs or
/// from two ends of an SV. These are returned as a vector of Connections.
/// Also gets locations of breakends (in any SVTYPE) from a VCF.
/// Stores them by variant ID in a hashmap of variant_id -> Vec<Coordinates>.
///
/// Returns:
///     vcf_breaks: vector of all unique breakend coordinates
///     coordinate_map: maps variant ID to a list of the genomic positions associated with it
pub fn get_vcf_breaks(
    vcf_filename: PathBuf,
    exclude_regions: &ExcludeRegions,
    target_region_opt: &Option<TargetCoordinate>,
) -> (Vec<Connection>, HashMap<String, Vec<Coordinate>>) {
    if should_skip_vcf_processing(&vcf_filename) {
        return (Vec::new(), HashMap::new());
    }

    let (mut bcf_reader, header) = open_vcf_reader(&vcf_filename);
    let mut result =
        process_vcf_records(&mut bcf_reader, &header, exclude_regions, target_region_opt);

    finalize_vcf_processing(&mut result, &header);

    (result.vcf_breaks, result.coordinate_map)
}

/// Checks if VCF processing should be skipped based on filename
fn should_skip_vcf_processing(vcf_filename: &Path) -> bool {
    let vcf_name_str = vcf_filename.as_os_str().to_str().unwrap();
    vcf_name_str == "\"\""
}

/// Opens VCF file and creates BCF reader with error handling
fn open_vcf_reader(vcf_filename: &PathBuf) -> (bcf::Reader, bcf::header::HeaderView) {
    let vcf_name_str = vcf_filename.as_os_str().to_str().unwrap();

    if !vcf_filename.exists() {
        error!("Error reading VCF {vcf_name_str}: file does not exist");
        std::process::exit(exitcode::IOERR);
    }

    let bcf_reader = match bcf::Reader::from_path(vcf_filename) {
        Ok(bcf) => bcf,
        Err(e) => {
            error!("Error reading {vcf_name_str}: {e}");
            std::process::exit(exitcode::IOERR);
        }
    };

    let header = bcf::Read::header(&bcf_reader).clone();
    (bcf_reader, header)
}

/// Processes all VCF records and populates the result structure
fn process_vcf_records(
    bcf_reader: &mut bcf::Reader,
    header: &bcf::header::HeaderView,
    exclude_regions: &ExcludeRegions,
    target_region_opt: &Option<TargetCoordinate>,
) -> VcfProcessingResult {
    let mut result = VcfProcessingResult::new();

    for record_result in bcf::Read::records(bcf_reader) {
        let record = record_result.expect("Fail to read record");
        let record_id = extract_record_id(&record);

        if should_skip_record(&record, header, exclude_regions, target_region_opt) {
            continue;
        }

        let svtype = get_optional_string_info(&record, "SVTYPE").unwrap_or_default();

        if svtype == "BND" {
            result.bnd_records.insert(record_id, record);
        } else if svtype != "INS" {
            process_standard_variant(&record, &record_id, header, &mut result);
        }
    }

    result
}

/// Extracts record ID from VCF record with error handling
fn extract_record_id(record: &bcf::Record) -> String {
    match String::from_utf8(record.id()) {
        Ok(id) => id,
        Err(e) => {
            error!("{e}");
            std::process::exit(exitcode::IOERR);
        }
    }
}

/// Determines if a record should be skipped based on exclusion criteria
fn should_skip_record(
    record: &bcf::Record,
    header: &bcf::header::HeaderView,
    exclude_regions: &ExcludeRegions,
    target_region_opt: &Option<TargetCoordinate>,
) -> bool {
    let chrom = get_vcf_record_chrom(record, header);
    let start = record.pos();
    let end = record.end();

    utils::entry_excluded(exclude_regions, (chrom, start, end), target_region_opt)
}

/// Processes standard (non-BND) variants and adds them to the result
fn process_standard_variant(
    record: &bcf::Record,
    record_id: &str,
    header: &bcf::header::HeaderView,
    result: &mut VcfProcessingResult,
) {
    let chrom = get_vcf_record_chrom(record, header);
    let start = record.pos();
    let end = record.end();

    let (pos_confidence_interval, end_confidence_interval) = get_confidence_intervals(record);
    let first_coord =
        Coordinate::new_with_confidence_interval(chrom.clone(), start, pos_confidence_interval);
    let second_coord =
        Coordinate::new_with_confidence_interval(chrom, end, end_confidence_interval);

    // Add first coordinate to mapping
    result
        .coordinate_map
        .entry(record_id.to_string())
        .or_default()
        .push(first_coord.clone());

    let dist = first_coord.start - second_coord.start;
    if dist.abs() > 1 {
        result
            .coordinate_map
            .entry(record_id.to_string())
            .or_default()
            .push(second_coord.clone());

        let connection = Connection::new(first_coord, second_coord, false, true);
        result.vcf_breaks.push(connection);
    }
}

/// Finalizes VCF processing by handling BND records and cleaning up coordinate map
fn finalize_vcf_processing(result: &mut VcfProcessingResult, header: &bcf::header::HeaderView) {
    add_bnd_connections(
        &mut result.vcf_breaks,
        &mut result.coordinate_map,
        std::mem::take(&mut result.bnd_records),
        header,
    );

    // Clean up coordinate map
    for (_, coords) in result.coordinate_map.iter_mut() {
        coords.sort();
        coords.dedup();
    }

    debug!(
        "{} unique variant coordinates found in VCF",
        result.vcf_breaks.len()
    );
}

/// Add the record that are BNDs in the vcf,
/// which are separate records connected by MATEID/ID
/// These are stored in a hashmap by ID -> record.
fn add_bnd_connections(
    vcf_breaks: &mut Vec<Connection>,
    coordinate_map: &mut HashMap<String, Vec<Coordinate>>,
    bnd_records: HashMap<String, bcf::Record>,
    header: &bcf::header::HeaderView,
) {
    for (first_bnd_id, first_bnd_record) in bnd_records.iter() {
        if let Some(mate_id) = get_optional_string_info(first_bnd_record, "MATEID") {
            if let Some(second_bnd_record) = bnd_records.get(&mate_id) {
                let first_chrom = get_vcf_record_chrom(first_bnd_record, header);
                let first_start = first_bnd_record.pos() + 1;
                let first_confidence_intervals = get_confidence_intervals(first_bnd_record);
                let first_record_coord = Coordinate::new_with_confidence_interval(
                    first_chrom.clone(),
                    first_start,
                    first_confidence_intervals.0,
                );

                let second_chrom = get_vcf_record_chrom(second_bnd_record, header);
                let second_start = second_bnd_record.pos() + 1;
                let second_confidence_intervals = get_confidence_intervals(first_bnd_record);
                let second_record_coord = Coordinate::new_with_confidence_interval(
                    second_chrom.clone(),
                    second_start,
                    second_confidence_intervals.0,
                );
                coordinate_map
                    .entry(first_bnd_id.clone())
                    .or_default()
                    .push(first_record_coord.clone());
                coordinate_map
                    .entry(mate_id.clone())
                    .or_default()
                    .push(second_record_coord.clone());

                vcf_breaks.push(Connection::new(
                    first_record_coord,
                    second_record_coord,
                    false,
                    true,
                ));
            }
        }
    }
}

/// Gets the chromosome name from a vcf record or dies trying.
/// Forces it to match the genome's expected chromosome if
/// a known human reference is chosen.
fn get_vcf_record_chrom(
    record: &rust_htslib::bcf::Record,
    header: &bcf::header::HeaderView,
) -> String {
    if let Some(ref_id) = record.rid() {
        if let Ok(chrom_bytes) = header.rid2name(ref_id) {
            let raw_chrom = String::from_utf8_lossy(chrom_bytes).to_string();
            let has_prefix = raw_chrom.get(0..3) == Some("chr");
            if has_prefix {
                return raw_chrom;
            } else {
                return format!("chr{raw_chrom}");
            }
        }
        error!(
            "Error reading VCF: Chromosome {ref_id} missing from header"
        );
        std::process::exit(exitcode::IOERR);
    }

    error!("Error reading VCF: Chromosome missing from header");
    std::process::exit(exitcode::IOERR);
}

/// Returns an String INFO field value from a VCF record
/// ## Arguments
/// * `record` - the variant record to check
/// * `field` - the name of the INFO field tag to return
fn get_optional_string_info(record: &rust_htslib::bcf::Record, field: &str) -> Option<String> {
    // check if this has the tag field and parse into an SV type if it does
    let mut tag_option_str = None;
    let tag_result = record.info(field.as_bytes()).string();
    if let Ok(Some(tag_value_vec)) = tag_result {
        // the tag_value_vec is an array of strings at this point, make sure we only get one
        if tag_value_vec.len() == 1 {
            if let Ok(tag_value_str) = String::from_utf8(tag_value_vec[0].to_vec()) {
                tag_option_str = Some(tag_value_str);
            }
        }
    }
    tag_option_str
}
/// Creates and return a tuple of positive confidence intervance from pos and end
/// using CIPOS and CIEND if present and greater than MAX_CLUST_DIST, else MAX_CLUST_DIST.
/// ## Arguments
/// * `record` - the variant record to check
fn get_confidence_intervals(record: &rust_htslib::bcf::Record) -> ((i64, i64), (i64, i64)) {
    // check if this has a CIPOS and CIEND fields and parse into integers if so
    let mut confidence_intervals = (
        (utils::MAX_CLUST_DIST, utils::MAX_CLUST_DIST),
        (utils::MAX_CLUST_DIST, utils::MAX_CLUST_DIST),
    );
    let cipos_result = record.info("CIPOS".as_bytes()).integer();
    if let Ok(Some(cipos)) = cipos_result {
        confidence_intervals.0 .0 = max(confidence_intervals.0 .0, cipos[0] as i64).abs();
        confidence_intervals.0 .1 = max(confidence_intervals.0 .1, cipos[1] as i64).abs();
    };

    let ciend_result = record.info("CIEND".as_bytes()).integer();
    if let Ok(Some(ciend)) = ciend_result {
        confidence_intervals.1 .0 = max(confidence_intervals.1 .0, ciend[0] as i64).abs();
        confidence_intervals.1 .1 = max(confidence_intervals.1 .1, ciend[1] as i64).abs();
    };
    confidence_intervals
}

/// If alignments are in reverse order, flip them.
/// Simplifies later processing conceptually.
/// If they aren't, return them unchanged.
fn flip_reverse_alignments(alignments: &mut [FwdStrandSplitReadSegment]) {
    let last_alignment = alignments.last().unwrap();
    // if the end of the last alignment is clipped,
    // this can't be the end of the alignment block
    if last_alignment.is_end_softclipped {
        alignments.reverse();
    }
}

/// Read the Sawfish json file that associates sample IDs and
/// readnames with variant IDs. Creates a {dict of sample ID ->{dict of variant ID -> vector of readnames}},
/// then returns only the dictionary for the given sample.
pub fn get_read_info_from_json(
    json_path: PathBuf,
    sample_id: String,
    variant_ids: HashSet<String>,
) -> HashMap<String, Vec<String>> {
    let json_name_str = json_path.as_os_str().to_str().unwrap();
    let file_handle = match std::fs::File::open(&json_path) {
        Ok(file) => file,
        Err(e) => {
            error!("Input JSON does not exist: \"{json_name_str}\"");
            error!("{e}");
            std::process::exit(exitcode::IOERR);
        }
    };
    let reader = BufReader::new(file_handle);
    let is_gzipped = is_gzipped(&json_name_str.to_string());

    let json_data_result: Result<HashMap<String, HashMap<String, Vec<String>>>, serde_json::Error> =
        if is_gzipped {
            let decoder = GzDecoder::new(reader);
            serde_json::from_reader(decoder)
        } else {
            serde_json::from_reader(reader)
        };

    let json_data: HashMap<String, HashMap<String, Vec<String>>> = match json_data_result {
        Ok(json_data) => json_data,
        Err(e) => {
            error!("Error reading JSON file \"{json_name_str}\"");
            error!("{e}");
            std::process::exit(exitcode::IOERR);
        }
    };
    let mut sample_variant_readnames = HashMap::new();
    for (variant_id, variant_data) in json_data.into_iter() {
        if variant_ids.contains(&variant_id) {
            if let Some(sample_json) = variant_data.get(&sample_id) {
                sample_variant_readnames.insert(variant_id, sample_json.to_owned());
            }
        }
    }
    if sample_variant_readnames.is_empty() {
        error!(
            "No variant data found for sample \"{sample_id}\" in \"{json_name_str}\""
        );
        std::process::exit(exitcode::DATAERR);
    }
    debug!(
        "{} variant IDs found in JSON",
        sample_variant_readnames.len()
    );
    sample_variant_readnames
}

/// Given a bedfile entry where a string is tab-delimited into
/// at least three fields, with the first being a string
/// and the next two being integers, create a coordinate to
/// represent them
fn bed_entry_to_coord(entry: String) -> Coordinate {
    let fields: Vec<&str> = entry.split_whitespace().collect();
    if fields.len() < 3 {
        error!("Invalid entry length in exclude regions file");
        std::process::exit(exitcode::IOERR);
    }
    let chrom = fields[0].to_string();

    let pos: i64 = match fields[1].parse() {
        Ok(num) => num,
        Err(_) => {
            error!("Non-numeric start position in exclude regions file");
            std::process::exit(exitcode::IOERR);
        }
    };
    let end: i64 = match fields[2].parse() {
        Ok(num) => num,
        Err(_) => {
            error!("Non-numeric end position in exclude regions file");
            std::process::exit(exitcode::IOERR);
        }
    };
    Coordinate::new_region(chrom, pos, end)
}

/// Reads a BED file of genomic regios to exclude. Can read a GZIP compressed or
/// uncompressed BAM from a local path. Regions are stored in an ExcludeRegions
/// struct with a HashMap of String chromosome name to Coordinate.
pub fn load_exclude_regions(exclude_regions_path: String) -> ExcludeRegions {
    let lines = if is_local_file(&exclude_regions_path) {
        utils::read_file_from_path(&exclude_regions_path)
    } else {
        error!(
            "Failed to read exclude regions file {exclude_regions_path}"
        );
        std::process::exit(exitcode::IOERR);
    };
    let mut exclude_regions = ExcludeRegions::new();
    for line in lines {
        if line.starts_with('#') {
            // skip headers
            continue;
        }
        let entry_coord = bed_entry_to_coord(line);
        exclude_regions
            .regions
            .entry(entry_coord.start_chrom.clone())
            .or_default()
            .insert(entry_coord);
    }
    exclude_regions
}
