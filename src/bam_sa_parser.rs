// from Chris Saunders, sawfish sv caller. Modified for complex sv

use log::error;
use rust_htslib::bam::{
    self,
    record::{Cigar, CigarString, Record},
};

use crate::{
    containers::{
        Coordinate, ExcludeRegions, FwdStrandSplitReadSegment, ReadMetadata,
        SegmentProcessingContext, SplitAlignmentSegment,
    },
    utils,
};
static MIN_MAPQ: u8 = 5;

fn parse_sa_segment(seg: &str) -> SplitAlignmentSegment {
    let sa_fields = seg.split_terminator(',').collect::<Vec<_>>();
    if sa_fields.len() != 6 {
        error!("Unexpected segment in bam SA tag: {seg}");
        std::process::exit(exitcode::IOERR);
    }
    let rname = sa_fields[0].to_string();
    let pos = sa_fields[1].parse::<i64>().unwrap() - 1;
    let is_fwd_strand = sa_fields[2] == "+";
    let cigar = CigarString::try_from(sa_fields[3].as_bytes()).unwrap();
    let mapq: u8 = sa_fields[4].parse::<u8>().unwrap();
    SplitAlignmentSegment {
        rname,
        pos,
        is_fwd_strand,
        cigar,
        mapq,
    }
}

fn parse_sa_aux_val(sa_aux_val: &str) -> Vec<SplitAlignmentSegment> {
    sa_aux_val
        .split_terminator(';')
        .map(parse_sa_segment)
        .collect::<Vec<SplitAlignmentSegment>>()
}

pub fn get_fwd_read_split_segments(
    record: &Record,
    chrom: String,
    exclude_regions: &mut ExcludeRegions,
) -> Vec<FwdStrandSplitReadSegment> {
    const SA_AUX_TAG: &[u8] = b"SA";
    const PS_AUX_TAG: &[u8] = b"PS";
    const HP_AUX_TAG: &[u8] = b"HP";
    let mut fwd_read_split_segments = Vec::new();
    let mut excluded = false;

    let readname = match String::from_utf8(record.qname().to_vec()) {
        Ok(name) => name,
        Err(_) => {
            error!("Failed to get read name from record");
            std::process::exit(exitcode::IOERR);
        }
    };

    // if mapq doesn't meet min threshold, skip the whole record
    let primary_mapq = record.mapq();
    if primary_mapq < MIN_MAPQ {
        return fwd_read_split_segments;
    }

    if let Some(sa_aux_val) = get_optional_string_aux_tag(record, SA_AUX_TAG) {
        let sa_segments = parse_sa_aux_val(&sa_aux_val);
        let phaseset_tag = get_optional_int_aux_tag(record, PS_AUX_TAG);
        let haplotype_tag = get_optional_int_aux_tag(record, HP_AUX_TAG);

        // Add the primary alignment first:
        let read_metadata = ReadMetadata {
            readname,
            phaseset_tag,
            haplotype_tag,
        };

        let mut context = SegmentProcessingContext {
            fwd_read_split_segments: &mut fwd_read_split_segments,
            exclude_regions,
            excluded: &mut excluded,
        };

        let primary_read_size =
            process_primary_alignment(record, &chrom, &read_metadata, &mut context);

        // Put the split and primary alignment segments in read
        // coordinates and sample order
        process_sa_segments(
            &sa_segments,
            primary_read_size,
            &read_metadata,
            &mut context,
        );
    }
    // if this has been excluded, add all alignments of the read to the exclude list
    if excluded {
        for alignment in fwd_read_split_segments.iter() {
            let exclude_coord =
                Coordinate::new_region(alignment.chrom.clone(), alignment.pos, alignment.end);
            exclude_regions
                .regions
                .entry(alignment.chrom.clone())
                .or_default()
                .insert(exclude_coord);
        }
        fwd_read_split_segments = Vec::new();
    }

    fwd_read_split_segments.sort_by_key(|x| x.fwd_read_start);
    generate_unspanned_segments(fwd_read_split_segments)
}

/// Generate alignment segments for reference connections between breaks
/// These don't have to be compared against the exclude regions because all of the
/// breakpoints are already checked when the spanned alignments are added.
fn generate_unspanned_segments(
    splits: Vec<FwdStrandSplitReadSegment>,
) -> Vec<FwdStrandSplitReadSegment> {
    let mut expanded_splits = Vec::new();
    let mut prev_chrom_opt: Option<String> = None;
    let mut prev_end_opt: Option<i64> = None;

    for split in splits {
        if let (Some(prev_chrom), Some(prev_end)) = (prev_chrom_opt, prev_end_opt) {
            let is_fwd_strand = true;
            let pos = prev_end;
            let mut end = split.pos;
            if !split.is_fwd_strand {
                end = split.end;
            }

            if prev_chrom < split.chrom {
                let unspanned_split = FwdStrandSplitReadSegment {
                    fwd_read_start: split.fwd_read_start,
                    fwd_read_end: split.fwd_read_start,
                    chrom: prev_chrom.clone(),
                    second_chrom: split.chrom.clone(),
                    pos,
                    end,
                    is_fwd_strand,
                    is_start_softclipped: true,
                    is_end_softclipped: true,
                    phaseset_tag: split.phaseset_tag,
                    haplotype_tag: split.haplotype_tag,
                    from_primary_bam_record: false,
                    readname: split.readname.clone(),
                    spans: false,
                };
                expanded_splits.push(unspanned_split);
            } else {
                let unspanned_split = FwdStrandSplitReadSegment {
                    fwd_read_start: split.fwd_read_start,
                    fwd_read_end: split.fwd_read_start,
                    chrom: split.chrom.clone(),
                    second_chrom: prev_chrom.clone(),
                    pos: end,
                    end: pos,
                    is_fwd_strand,
                    is_start_softclipped: true,
                    is_end_softclipped: true,
                    phaseset_tag: split.phaseset_tag,
                    haplotype_tag: split.haplotype_tag,
                    from_primary_bam_record: false,
                    readname: split.readname.clone(),
                    spans: false,
                };
                expanded_splits.push(unspanned_split);
            }
        }

        prev_chrom_opt = Some(split.chrom.clone());
        prev_end_opt = Some(split.end);
        if !split.is_fwd_strand {
            prev_end_opt = Some(split.pos);
        }
        expanded_splits.push(split);
    }
    expanded_splits
}

fn get_sa_softclip_status(cigar: &CigarString, min_softclips: i64) -> (bool, bool) {
    let is_start_softclipped = match cigar.first().unwrap() {
        Cigar::SoftClip(clip_len) => *clip_len as i64 > min_softclips,
        _ => false,
    };
    let is_end_softclipped = match cigar.last().unwrap() {
        Cigar::SoftClip(clip_len) => *clip_len as i64 > min_softclips,
        _ => false,
    };
    (is_start_softclipped, is_end_softclipped)
}

fn get_fwd_strand_read_pos(
    read_start: usize,
    read_end: usize,
    read_size: usize,
    is_fwd_strand: bool,
) -> (usize, usize) {
    if is_fwd_strand {
        (read_start, read_end)
    } else {
        (read_size - read_end, read_size - read_start)
    }
}

/// Retrieve a string aux tag from bam file
///
/// Function will panic if the tag has a non-string value
///
pub fn get_optional_string_aux_tag(record: &bam::Record, aux_tag: &[u8]) -> Option<String> {
    match record.aux(aux_tag) {
        Ok(aux_val) => Some(match aux_val {
            bam::record::Aux::String(val) => val.to_string(),
            _ => unexpected_aux_val_err(record, aux_tag, aux_val),
        }),
        _ => None,
    }
}

/// Retrieve an integer aux tag from bam file
///
/// Function will panic if the tag has a non-int value
///
pub fn get_optional_int_aux_tag(record: &bam::Record, aux_tag: &[u8]) -> Option<i32> {
    match record.aux(aux_tag) {
        Ok(aux_val) => Some(match aux_val {
            bam::record::Aux::U32(val) => val as i32,
            bam::record::Aux::I32(val) => val,
            bam::record::Aux::U16(val) => val as i32,
            bam::record::Aux::I16(val) => val as i32,
            bam::record::Aux::U8(val) => val as i32,
            bam::record::Aux::I8(val) => val as i32,
            _ => unexpected_aux_val_err(record, aux_tag, aux_val),
        }),
        _ => None,
    }
}

/// Report the following positions in read coordinates:
/// 1. The first position after all left-side hard and soft clipping
/// 2. The first position before all right-side hard and soft clipping
/// 3. The read length
///
fn get_complete_read_clip_positions(cigar: &[Cigar]) -> (usize, usize, usize) {
    use Cigar::*;

    let mut ref_pos = 0;
    let mut read_pos = 0;
    let mut left_clip_size = 0;
    let mut right_clip_size = 0;
    let mut left_clip = true;
    for c in cigar.iter() {
        match c {
            HardClip(len) | SoftClip(len) => {
                if left_clip {
                    left_clip_size += *len as usize;
                } else {
                    right_clip_size += *len as usize;
                }
            }
            _ => {
                left_clip = false;
            }
        };
        update_ref_and_complete_read_pos(c, &mut ref_pos, &mut read_pos);
    }
    (left_clip_size, read_pos - right_clip_size, read_pos)
}

/// A utility method to track ref and read positions while iterating through a cigar string
///
/// Read position here is the position in the original read, before hard-clipping
///
/// # Example
/// ```ignore
/// let mut ref_pos = 100;
/// let mut read_pos = 100;
/// for (index, c) in record.cigar().iter().enumerate() {
///     update_ref_and_complete_read_pos(c, &mut ref_pos, &mut read_pos);
/// }
/// ```
fn update_ref_and_complete_read_pos(c: &Cigar, ref_pos: &mut i64, read_pos: &mut usize) {
    update_complete_read_pos(c, read_pos);
    update_ref_pos(c, ref_pos);
}

/// A utility method to track read positions while iterating through a cigar string
///
/// Read position here is the position in the original read, before hard-clipping
fn update_complete_read_pos(c: &Cigar, read_pos: &mut usize) {
    *read_pos += get_cigarseg_complete_read_offset(c);
}

/// A utility method to track ref positions while iterating through a cigar string
fn update_ref_pos(c: &Cigar, ref_pos: &mut i64) {
    *ref_pos += get_cigarseg_ref_offset(c);
}

fn get_cigarseg_ref_offset(c: &Cigar) -> i64 {
    use Cigar::*;
    match c {
        Del(len) | RefSkip(len) | Diff(len) | Equal(len) | Match(len) => *len as i64,
        _ => 0,
    }
}

fn get_cigarseg_complete_read_offset(c: &Cigar) -> usize {
    use Cigar::*;
    match c {
        HardClip(len) | Ins(len) | SoftClip(len) | Diff(len) | Equal(len) | Match(len) => {
            *len as usize
        }
        _ => 0,
    }
}

fn unexpected_aux_val_err(
    record: &bam::Record,
    aux_tag: &[u8],
    aux_val: bam::record::Aux<'_>,
) -> ! {
    let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
    panic!(
        "Unexpected {} tag format in read {qname}: {:?}",
        std::str::from_utf8(aux_tag).unwrap(),
        aux_val,
    );
}

fn get_seq_len_from_cigar(cigar: &CigarString) -> i64 {
    let mut aligned_len = 0;
    for operation in cigar.iter() {
        match operation {
            Cigar::Match(len)
            | Cigar::Equal(len)
            | Cigar::Del(len)
            | Cigar::RefSkip(len)
            | Cigar::Diff(len) => aligned_len += len,
            _ => {}
        }
    }
    aligned_len as i64
}

/// Process the primary alignment record and create a FwdStrandSplitReadSegment.
/// Checks for exclusion regions and returns the total read size for validation.
fn process_primary_alignment(
    record: &Record,
    chrom: &str,
    read_metadata: &ReadMetadata,
    context: &mut SegmentProcessingContext,
) -> usize {
    const MIN_SOFTCLIPS: i64 = 10;
    let (read_start, read_end, read_size) = get_complete_read_clip_positions(&record.cigar());
    let (fwd_read_start, fwd_read_end) =
        get_fwd_strand_read_pos(read_start, read_end, read_size, !record.is_reverse());

    let aligned_len = get_seq_len_from_cigar(&record.cigar().take());
    let is_start_softclipped = record.cigar().leading_softclips() > MIN_SOFTCLIPS;
    let is_end_softclipped = record.cigar().trailing_softclips() > MIN_SOFTCLIPS;
    let alignment = FwdStrandSplitReadSegment {
        fwd_read_start,
        fwd_read_end,
        chrom: chrom.to_string(),
        second_chrom: chrom.to_string(),
        pos: record.pos(),
        end: record.pos() + aligned_len,
        is_fwd_strand: !record.is_reverse(),
        is_start_softclipped,
        is_end_softclipped,
        phaseset_tag: read_metadata.phaseset_tag,
        haplotype_tag: read_metadata.haplotype_tag,
        from_primary_bam_record: true,
        readname: read_metadata.readname.clone(),
        spans: true,
    };
    let entry = (alignment.chrom.clone(), alignment.pos, alignment.end);
    if !*context.excluded {
        *context.excluded = utils::entry_excluded(context.exclude_regions, entry, &None);
    }
    context.fwd_read_split_segments.push(alignment);
    read_size
}

/// Process supplementary alignment segments from the SA tag.
/// Creates FwdStrandSplitReadSegment objects for each valid SA segment and handles MAPQ filtering.
fn process_sa_segments(
    sa_segments: &[SplitAlignmentSegment],
    primary_read_size: usize,
    read_metadata: &ReadMetadata,
    context: &mut SegmentProcessingContext,
) {
    const MIN_SOFTCLIPS: i64 = 10;
    for sa_segment in sa_segments.iter() {
        if sa_segment.mapq < MIN_MAPQ {
            *context.excluded = true;
        }
        let (read_start, read_end, read_size) = get_complete_read_clip_positions(&sa_segment.cigar);
        assert_eq!(primary_read_size, read_size);
        let (fwd_read_start, fwd_read_end) =
            get_fwd_strand_read_pos(read_start, read_end, read_size, sa_segment.is_fwd_strand);
        let aligned_len = get_seq_len_from_cigar(&sa_segment.cigar);
        if sa_segment.cigar.is_empty() {
            continue;
        }
        let (is_start_softclipped, is_end_softclipped) =
            get_sa_softclip_status(&sa_segment.cigar, MIN_SOFTCLIPS);
        let alignment = FwdStrandSplitReadSegment {
            fwd_read_start,
            fwd_read_end,
            chrom: sa_segment.rname.clone(),
            second_chrom: sa_segment.rname.clone(),
            phaseset_tag: read_metadata.phaseset_tag,
            haplotype_tag: read_metadata.haplotype_tag,
            pos: sa_segment.pos,
            end: sa_segment.pos + aligned_len,
            is_fwd_strand: sa_segment.is_fwd_strand,
            is_start_softclipped,
            is_end_softclipped,
            from_primary_bam_record: false,
            readname: read_metadata.readname.clone(),
            spans: true,
        };
        let entry = (alignment.chrom.clone(), alignment.pos, alignment.end);
        if !*context.excluded {
            *context.excluded = utils::entry_excluded(context.exclude_regions, entry, &None);
        }
        context.fwd_read_split_segments.push(alignment);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;

    // Helper function to create a test CIGAR string
    fn create_test_cigar(cigar_str: &str) -> CigarString {
        CigarString::try_from(cigar_str.as_bytes()).unwrap()
    }

    #[test]
    fn test_parse_sa_segment() {
        let test_segment = "chr3,10001,+,5535S10=1D39=2X11438S,60,192";
        let result = parse_sa_segment(test_segment);

        assert_eq!(result.rname, "chr3");
        assert_eq!(result.pos, 10000); // Should be 0-indexed
        assert!(result.is_fwd_strand);
        assert_eq!(result.mapq, 60);
        assert_eq!(result.cigar.to_string(), "5535S10=1D39=2X11438S");
    }

    #[test]
    fn test_parse_sa_segment_reverse_strand() {
        let test_segment = "chr4,106872270,-,23=1I226=1I195=1X147=1D1021=7362S,60,19";
        let result = parse_sa_segment(test_segment);

        assert_eq!(result.rname, "chr4");
        assert_eq!(result.pos, 106872269); // Should be 0-indexed
        assert!(!result.is_fwd_strand);
        assert_eq!(result.mapq, 60);
    }

    // Note: We can't easily test the invalid format case since parse_sa_segment
    // calls std::process::exit instead of panicking

    #[test]
    fn test_parse_sa_aux_val() {
        let test_val = "chr3,10001,+,5535S10=1D39=2X11438S,60,192;\
        chr3,10001,+,3073S15=2D20=2X11=1X5=1I23=1X5=14798S,22,44;\
        chr4,106872270,-,23=1I226=1I195=1X147=1D1021=7362S,60,19;";

        let result = parse_sa_aux_val(test_val);

        assert_eq!(result.len(), 3);
        assert_eq!(result[2].rname, "chr4");
        assert_eq!(result[1].pos, 10_000);
        assert!(!result[2].is_fwd_strand);
    }

    #[test]
    fn test_parse_sa_aux_val_empty() {
        let result = parse_sa_aux_val("");
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_parse_sa_aux_val_single_segment() {
        let test_val = "chr1,1000,+,100M,30,1;";
        let result = parse_sa_aux_val(test_val);

        assert_eq!(result.len(), 1);
        assert_eq!(result[0].rname, "chr1");
        assert_eq!(result[0].pos, 999); // 0-indexed
    }

    #[test]
    fn test_get_sa_softclip_status() {
        // Test both ends softclipped with sufficient length
        let cigar = create_test_cigar("15S100M20S");
        let (start_clipped, end_clipped) = get_sa_softclip_status(&cigar, 10);
        assert!(start_clipped);
        assert!(end_clipped);

        // Test only start softclipped
        let cigar = create_test_cigar("15S100M");
        let (start_clipped, end_clipped) = get_sa_softclip_status(&cigar, 10);
        assert!(start_clipped);
        assert!(!end_clipped);

        // Test insufficient softclips
        let cigar = create_test_cigar("5S100M5S");
        let (start_clipped, end_clipped) = get_sa_softclip_status(&cigar, 10);
        assert!(!start_clipped);
        assert!(!end_clipped);

        // Test no softclips
        let cigar = create_test_cigar("100M");
        let (start_clipped, end_clipped) = get_sa_softclip_status(&cigar, 10);
        assert!(!start_clipped);
        assert!(!end_clipped);
    }

    #[test]
    fn test_get_fwd_strand_read_pos_forward() {
        let (fwd_start, fwd_end) = get_fwd_strand_read_pos(10, 90, 100, true);
        assert_eq!(fwd_start, 10);
        assert_eq!(fwd_end, 90);
    }

    #[test]
    fn test_get_fwd_strand_read_pos_reverse() {
        let (fwd_start, fwd_end) = get_fwd_strand_read_pos(10, 90, 100, false);
        assert_eq!(fwd_start, 10); // 100 - 90 = 10
        assert_eq!(fwd_end, 90); // 100 - 10 = 90
    }

    #[test]
    fn test_get_complete_read_clip_positions() {
        use Cigar::*;

        // Test with soft and hard clips on both ends
        let cigar = vec![SoftClip(10), HardClip(5), Match(80), SoftClip(15)];
        let (left_clip, right_clip, total) = get_complete_read_clip_positions(&cigar);
        assert_eq!(left_clip, 15); // 10 + 5
        assert_eq!(right_clip, 95); // 110 - 15
        assert_eq!(total, 110); // 10 + 5 + 80 + 15

        // Test with no clips
        let cigar = vec![Match(100)];
        let (left_clip, right_clip, total) = get_complete_read_clip_positions(&cigar);
        assert_eq!(left_clip, 0);
        assert_eq!(right_clip, 100);
        assert_eq!(total, 100);

        // Test with only left clips
        let cigar = vec![SoftClip(20), Match(80)];
        let (left_clip, right_clip, total) = get_complete_read_clip_positions(&cigar);
        assert_eq!(left_clip, 20);
        assert_eq!(right_clip, 100);
        assert_eq!(total, 100);
    }

    #[test]
    fn test_get_cigarseg_ref_offset() {
        use Cigar::*;

        assert_eq!(get_cigarseg_ref_offset(&Match(50)), 50);
        assert_eq!(get_cigarseg_ref_offset(&Equal(30)), 30);
        assert_eq!(get_cigarseg_ref_offset(&Diff(20)), 20);
        assert_eq!(get_cigarseg_ref_offset(&Del(10)), 10);
        assert_eq!(get_cigarseg_ref_offset(&RefSkip(100)), 100);

        // Operations that don't consume reference
        assert_eq!(get_cigarseg_ref_offset(&Ins(25)), 0);
        assert_eq!(get_cigarseg_ref_offset(&SoftClip(40)), 0);
        assert_eq!(get_cigarseg_ref_offset(&HardClip(15)), 0);
    }

    #[test]
    fn test_get_cigarseg_complete_read_offset() {
        use Cigar::*;

        assert_eq!(get_cigarseg_complete_read_offset(&Match(50)), 50);
        assert_eq!(get_cigarseg_complete_read_offset(&Equal(30)), 30);
        assert_eq!(get_cigarseg_complete_read_offset(&Diff(20)), 20);
        assert_eq!(get_cigarseg_complete_read_offset(&Ins(25)), 25);
        assert_eq!(get_cigarseg_complete_read_offset(&SoftClip(40)), 40);
        assert_eq!(get_cigarseg_complete_read_offset(&HardClip(15)), 15);

        // Operations that don't consume read
        assert_eq!(get_cigarseg_complete_read_offset(&Del(10)), 0);
        assert_eq!(get_cigarseg_complete_read_offset(&RefSkip(100)), 0);
    }

    #[test]
    fn test_update_ref_pos() {
        use Cigar::*;

        let mut ref_pos = 100i64;
        update_ref_pos(&Match(50), &mut ref_pos);
        assert_eq!(ref_pos, 150);

        update_ref_pos(&Del(10), &mut ref_pos);
        assert_eq!(ref_pos, 160);

        update_ref_pos(&Ins(20), &mut ref_pos);
        assert_eq!(ref_pos, 160); // No change for insertion
    }

    #[test]
    fn test_update_complete_read_pos() {
        use Cigar::*;

        let mut read_pos = 50usize;
        update_complete_read_pos(&Match(30), &mut read_pos);
        assert_eq!(read_pos, 80);

        update_complete_read_pos(&Ins(10), &mut read_pos);
        assert_eq!(read_pos, 90);

        update_complete_read_pos(&Del(15), &mut read_pos);
        assert_eq!(read_pos, 90); // No change for deletion
    }

    #[test]
    fn test_update_ref_and_complete_read_pos() {
        use Cigar::*;

        let mut ref_pos = 100i64;
        let mut read_pos = 50usize;

        update_ref_and_complete_read_pos(&Match(30), &mut ref_pos, &mut read_pos);
        assert_eq!(ref_pos, 130);
        assert_eq!(read_pos, 80);

        update_ref_and_complete_read_pos(&Ins(10), &mut ref_pos, &mut read_pos);
        assert_eq!(ref_pos, 130); // No change in ref
        assert_eq!(read_pos, 90); // Read advances

        update_ref_and_complete_read_pos(&Del(5), &mut ref_pos, &mut read_pos);
        assert_eq!(ref_pos, 135); // Ref advances
        assert_eq!(read_pos, 90); // No change in read
    }

    #[test]
    fn test_get_seq_len_from_cigar() {
        let cigar = create_test_cigar("10S50M5I10M20D30M10S");
        let result = get_seq_len_from_cigar(&cigar);
        assert_eq!(result, 110); // 50 + 10 + 20 + 30 = 110 (excluding soft clips and insertions)

        let simple_cigar = create_test_cigar("100M");
        let result = get_seq_len_from_cigar(&simple_cigar);
        assert_eq!(result, 100);

        let complex_cigar = create_test_cigar("20=10X5D15=");
        let result = get_seq_len_from_cigar(&complex_cigar);
        assert_eq!(result, 50); // 20 + 10 + 5 + 15 = 50
    }

    #[test]
    fn test_generate_unspanned_segments() {
        // Create test split segments
        let splits = vec![
            create_test_split_segment("chr1", 100, 200, 0, 100, true),
            create_test_split_segment("chr1", 300, 400, 100, 200, true),
            create_test_split_segment("chr2", 500, 600, 200, 300, true),
        ];

        let result = generate_unspanned_segments(splits);

        // Should have original 3 segments plus 2 unspanned connections
        assert_eq!(result.len(), 5);

        // Check that unspanned segments are created
        let unspanned_segments: Vec<_> = result.iter().filter(|s| !s.spans).collect();
        assert_eq!(unspanned_segments.len(), 2);

        // First unspanned connects chr1 segments
        // Since "chr1" == "chr1" (else branch), pos=end=300, end=pos=200
        assert_eq!(unspanned_segments[0].chrom, "chr1");
        assert_eq!(unspanned_segments[0].second_chrom, "chr1");
        assert_eq!(unspanned_segments[0].pos, 300);
        assert_eq!(unspanned_segments[0].end, 200);

        // Second unspanned connects chr1 to chr2
        // Since "chr1" < "chr2", uses first branch: pos = prev_end (400), end = split.pos (500)
        assert_eq!(unspanned_segments[1].chrom, "chr1");
        assert_eq!(unspanned_segments[1].second_chrom, "chr2");
        assert_eq!(unspanned_segments[1].pos, 400);
        assert_eq!(unspanned_segments[1].end, 500);
    }

    #[test]
    fn test_generate_unspanned_segments_empty() {
        let result = generate_unspanned_segments(vec![]);
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_generate_unspanned_segments_single() {
        let splits = vec![create_test_split_segment("chr1", 100, 200, 0, 100, true)];

        let result = generate_unspanned_segments(splits);
        assert_eq!(result.len(), 1); // Only the original segment
        assert!(result[0].spans);
    }

    #[test]
    fn test_generate_unspanned_segments_reverse_strand() {
        let splits = vec![
            create_test_split_segment("chr1", 100, 200, 0, 100, true),
            create_test_split_segment("chr1", 300, 400, 100, 200, false), // Reverse strand
        ];

        let result = generate_unspanned_segments(splits);
        assert_eq!(result.len(), 3); // 2 original + 1 unspanned

        let unspanned = result.iter().find(|s| !s.spans).unwrap();
        // Since "chr1" == "chr1" (else branch), pos=end=400 (split.end for reverse), end=pos=200 (prev_end)
        assert_eq!(unspanned.pos, 400);
        assert_eq!(unspanned.end, 200);
    }

    // Helper function to create test FwdStrandSplitReadSegment
    fn create_test_split_segment(
        chrom: &str,
        pos: i64,
        end: i64,
        fwd_start: usize,
        fwd_end: usize,
        is_fwd_strand: bool,
    ) -> FwdStrandSplitReadSegment {
        FwdStrandSplitReadSegment {
            fwd_read_start: fwd_start,
            fwd_read_end: fwd_end,
            chrom: chrom.to_string(),
            second_chrom: chrom.to_string(),
            pos,
            end,
            is_fwd_strand,
            is_start_softclipped: false,
            is_end_softclipped: false,
            phaseset_tag: None,
            haplotype_tag: None,
            from_primary_bam_record: true,
            readname: "test_read".to_string(),
            spans: true,
        }
    }

    // Additional tests for edge cases and integration scenarios
    #[test]
    fn test_cigar_position_tracking_complex() {
        use Cigar::*;

        // Test a complex CIGAR with multiple operation types
        let cigar = vec![
            SoftClip(10), // 0-9 (read pos advances)
            Match(20),    // 10-29 (both advance)
            Ins(5),       // 30-34 (read advances, ref stays)
            Del(8),       // read stays, ref advances 8
            Match(15),    // 35-49 (both advance)
            SoftClip(5),  // 50-54 (read advances)
        ];

        let mut ref_pos = 100i64;
        let mut read_pos = 0usize;

        for c in &cigar {
            update_ref_and_complete_read_pos(c, &mut ref_pos, &mut read_pos);
        }

        assert_eq!(ref_pos, 143); // 100 + 20 + 8 + 15 = 143
        assert_eq!(read_pos, 55); // 10 + 20 + 5 + 15 + 5 = 55
    }

    #[test]
    fn test_softclip_edge_cases() {
        // Test minimum threshold
        let cigar = create_test_cigar("10S100M10S");
        let (start_clipped, end_clipped) = get_sa_softclip_status(&cigar, 10);
        assert!(!start_clipped); // Exactly at threshold, not greater
        assert!(!end_clipped);

        // Test just above threshold
        let cigar = create_test_cigar("11S100M11S");
        let (start_clipped, end_clipped) = get_sa_softclip_status(&cigar, 10);
        assert!(start_clipped);
        assert!(end_clipped);
    }

    #[test]
    fn test_fwd_strand_position_edge_cases() {
        // Test when read positions are at extremes
        let (fwd_start, fwd_end) = get_fwd_strand_read_pos(0, 100, 100, true);
        assert_eq!(fwd_start, 0);
        assert_eq!(fwd_end, 100);

        let (fwd_start, fwd_end) = get_fwd_strand_read_pos(0, 100, 100, false);
        assert_eq!(fwd_start, 0); // 100 - 100 = 0
        assert_eq!(fwd_end, 100); // 100 - 0 = 100
    }

    #[test]
    fn test_seq_len_complex_cigar() {
        // Test CIGAR with all types of operations
        let cigar = create_test_cigar("5H10S20M5I15D10=5X20N25M10S5H");
        let result = get_seq_len_from_cigar(&cigar);
        // Should sum: 20M + 15D + 10= + 5X + 20N + 25M = 95
        assert_eq!(result, 95);
    }

    #[test]
    fn test_generate_unspanned_chromosomal_order() {
        // Test that chromosomes are ordered correctly (chr1 < chr2)
        let splits = vec![
            create_test_split_segment("chr2", 100, 200, 0, 100, true),
            create_test_split_segment("chr1", 300, 400, 100, 200, true),
        ];

        let result = generate_unspanned_segments(splits);
        let unspanned = result.iter().find(|s| !s.spans).unwrap();

        // Since "chr2" > "chr1", uses second branch: chrom = split.chrom (chr1), second_chrom = prev_chrom (chr2)
        // pos = end (300), end = pos (200)
        assert_eq!(unspanned.chrom, "chr1");
        assert_eq!(unspanned.second_chrom, "chr2");
        assert_eq!(unspanned.pos, 300);
        assert_eq!(unspanned.end, 200);
    }
}
