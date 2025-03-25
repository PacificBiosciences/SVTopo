// from Chris Saunders, sawfish sv caller. Modified for complex sv

use log::error;
use rust_htslib::bam::{
    self,
    record::{Cigar, CigarString, Record},
};

use crate::{
    containers::{Coordinate, ExcludeRegions, FwdStrandSplitReadSegment, SplitAlignmentSegment},
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
    const MIN_SOFTCLIPS: i64 = 10;
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
        let primary_read_size = {
            let (read_start, read_end, read_size) =
                get_complete_read_clip_positions(&record.cigar());
            let (fwd_read_start, fwd_read_end) =
                get_fwd_strand_read_pos(read_start, read_end, read_size, !record.is_reverse());

            let aligned_len = get_seq_len_from_cigar(&record.cigar().take());
            let is_start_softclipped = record.cigar().leading_softclips() > MIN_SOFTCLIPS;
            let is_end_softclipped = record.cigar().trailing_softclips() > MIN_SOFTCLIPS;
            let alignment = FwdStrandSplitReadSegment {
                fwd_read_start,
                fwd_read_end,
                chrom: chrom.clone(),
                second_chrom: chrom.clone(),
                pos: record.pos(),
                end: record.pos() + aligned_len,
                is_fwd_strand: !record.is_reverse(),
                is_start_softclipped,
                is_end_softclipped,
                phaseset_tag,
                haplotype_tag,
                from_primary_bam_record: true,
                readname: readname.clone(),
                spans: true,
            };
            let entry = (alignment.chrom.clone(), alignment.pos, alignment.end);
            if !excluded {
                excluded = utils::entry_excluded(exclude_regions, entry, &None);
            }
            fwd_read_split_segments.push(alignment);
            read_size
        };

        // Put the split and primary alignment segments in read
        // coordinates and sample order
        for sa_segment in sa_segments.iter() {
            // if mapq doesn't meet min threshold, skip the whole record
            if sa_segment.mapq < MIN_MAPQ {
                excluded = true;
            }
            let (read_start, read_end, read_size) =
                get_complete_read_clip_positions(&sa_segment.cigar);
            assert_eq!(primary_read_size, read_size);
            let (fwd_read_start, fwd_read_end) =
                get_fwd_strand_read_pos(read_start, read_end, read_size, sa_segment.is_fwd_strand);
            let aligned_len = get_seq_len_from_cigar(&sa_segment.cigar);
            if sa_segment.cigar.len() == 0 {
                continue;
            }
            let (is_start_softclipped, is_end_softclipped) =
                get_sa_softclip_status(&sa_segment.cigar, MIN_SOFTCLIPS);
            let alignment = FwdStrandSplitReadSegment {
                fwd_read_start,
                fwd_read_end,
                chrom: sa_segment.rname.clone(),
                second_chrom: sa_segment.rname.clone(),
                phaseset_tag,
                haplotype_tag,
                pos: sa_segment.pos,
                end: sa_segment.pos + aligned_len,
                is_fwd_strand: sa_segment.is_fwd_strand,
                is_start_softclipped,
                is_end_softclipped,
                from_primary_bam_record: false,
                readname: readname.clone(),
                spans: true,
            };
            let entry = (alignment.chrom.clone(), alignment.pos, alignment.end);
            if !excluded {
                excluded = utils::entry_excluded(exclude_regions, entry, &None);
            }
            fwd_read_split_segments.push(alignment);
        }
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

#[cfg(test)]
mod tests {
    use super::*;

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
}
