#!/usr/bin/env python

import re
import os
import sys
import json
import gzip
import logging
from collections import namedtuple

logger = logging.getLogger(__name__)

# Define all namedtuples at module level
Region = namedtuple("Region", ["chrom", "start", "end", "size"])
BedRecord = namedtuple("BedRecord", ["start", "end", "strand", "title"])
AnnotationRecord = namedtuple(
    "AnnotationRecord", ["chrom", "start", "end", "title", "strand"]
)

ARROWHEAD_SIZE = 50
WINDOW_EDGE_SIZE = 20
MIN_CLUSTER_PRECISION = 10
SMALL_BLOCK_SCALE = 30
FONTSIZE = 7
SMALL_FONTSIZE = 4
PATCH_HALFWIDTH = 0.25
IMAGE_HEIGHT = 5
IMAGE_WIDTH = 12
CHAINS = "chains"  # key for axis lookup
BREAKS = "breaks"  # key for axis lookup
MEGABASE = 1_000_000
SPLIT_REGION_PADDING = 1_000  # padding for region splitting
BIG_NUMBER = 1_000_000_000_000_000_000
MAX_XTICK_COUNT = 8
MAX_ALLOWED_COVERAGE = 15
COVERAGE_COLOR = "black"
UNIQUE_BOX_COLOR = "darkgrey"
AMBIG_BOX_COLOR = "red"


def order_by_sample_idx(region_info):
    """
    from a list of region info, creates a dictionary of info entries
    by their order in the sample. Return the dictionary and the max index
    """
    reordered_region_info = {}
    sample_order_max = 0
    for block in region_info:
        if block["sample_order_index"] > sample_order_max:
            sample_order_max = block["sample_order_index"]
        if block["sample_order_index"] not in reordered_region_info:
            reordered_region_info[block["sample_order_index"]] = []
        reordered_region_info[block["sample_order_index"]].append(block)
    return reordered_region_info, sample_order_max


def get_coordinates(region_window: str):
    """
    Converts a string with chr:pos-end to namedtuple
    """
    try:
        window_components = region_window.split(":")
        if len(window_components) != 2:
            raise ValueError(f"Malformed region string: {region_window}")
        chrom = window_components[0][0].upper() + window_components[0][1:]
        if "-" in window_components[1]:
            start, end = window_components[1].split("-")
        else:
            start, end = window_components[1], window_components[1]
        chrom = chrom.upper()
        start = int(start)
        end = int(end)
        if start > end:
            start, end = end, start
        return Region(chrom, start, end, end - start)
    except (ValueError, IndexError):
        sys.exit("failed to parse region " + region_window)


def natural_sort_key(string_val):
    """
    sorts strings by number if any found, and alphabetically if not
    """
    result = []
    for text in re.split("([0-9]+)", string_val):
        if text.isdigit():
            result.append(int(text))
        else:
            result.append(text.upper())
    return result


def get_full_region(region_info):
    """
    Gets the full genomic region or window from a list of region info entries
    """
    positions = []

    for entry in region_info:
        for coverage in entry["coverages"]:
            positions.append(coverage)
        start = "{}:{}".format(entry["region"]["start_chrom"], entry["region"]["start"])
        end = "{}:{}".format(entry["region"]["end_chrom"], entry["region"]["end"])
        positions.append(start)
        positions.append(end)
    positions = sorted(positions, key=natural_sort_key)

    if len(positions) == 0:
        return
    elif len(positions) == 1:
        return get_coordinates(positions[0])
    else:
        start = get_coordinates(positions[0])
        end = get_coordinates(positions[-1])
        # Create a new Region with the start coordinates from start and end coordinates from end
        return Region(start.chrom, start.start, end.end, end.end - start.start)


def round_genomic_coords_to_str(start, end, max_at_kb=False):
    try:
        start = int(start)
    except (ValueError, OverflowError) as e:
        logger.error("ERROR: Region start is non-integer")
        logger.error(e)
        sys.exit(1)
    try:
        end = int(end)
    except (ValueError, OverflowError) as e:
        logger.error("ERROR: Region end is non-integer")
        logger.error(e)
        sys.exit(1)

    length = abs(end - start)
    len_str = ""
    if length < 1_000:
        len_str = "{:,} bp".format(length)
    elif 1_000 <= length <= 10_000:
        # if greater than 1kb, round to nearest hundred and report as kb
        rounded_length = int(round(length / 1_00) / 10)
        len_str = "{:,} kbp".format(rounded_length)
    elif 10_000 < length < 1_000_000 or max_at_kb:
        # if greater than 10kb, round to nearest thousand and report as kb
        rounded_length = round(length / 1_000)
        len_str = "{:,} kbp".format(rounded_length)
    else:
        # if greater than 1mb, round to nearest million and report as mb
        rounded_length = round(length / 1_000_000)
        len_str = "{:,} mbp".format(rounded_length)
    return len_str


def is_gzipped(putative_zipfile):
    """
    Check if file is zipped
    """
    with open(putative_zipfile, "rb") as filehandle:
        id_bytes = filehandle.read(2)
        return id_bytes == b"\x1f\x8b"


def unpack_json(json_filename):
    if json_filename.endswith("gz"):
        json_prefix = os.path.splitext(
            os.path.splitext(os.path.basename(json_filename))[0]
        )[0]
    else:
        json_prefix = os.path.splitext(os.path.basename(json_filename))[0]

    if not (os.path.exists(json_filename) and os.path.isfile(json_filename)):
        logger.warning(" {} does not exist".format(json_filename))
        return
    if json_filename.endswith(".gz"):
        if not is_gzipped(json_filename):
            logger.warning(
                "{} is identified as gzipped but is not".format(json_filename)
            )
            return
        with gzip.open(json_filename, "rt", encoding="UTF-8") as json_fh:
            try:
                return json.load(json_fh), json_prefix
            except json.decoder.JSONDecodeError:
                logger.warning(" {} is empty or misformatted".format(json_filename))
                return
    else:
        with open(json_filename, "r") as json_fh:
            try:
                return json.load(json_fh), json_prefix
            except json.decoder.JSONDecodeError:
                logger.warning(" {} is empty or misformatted".format(json_filename))
                return


def get_annotation_from_bed_record(bed_record, filename, line_number, _bed):
    """
    Get the necessary info from a line extracted from a bed file,
    including strand and record title if present. Return as
    annotation_record.
    """
    fields = bed_record.strip().split()
    if len(fields) < 3:
        logger.error(
            "Annotation file {} contains invalid record at line {}".format(
                filename, line_number
            )
        )
        sys.exit(1)
    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    # enforce start less than end
    if start > end:
        tmp = start
        start = end
        end = tmp
    title = ""
    if len(fields) > 3:
        title = fields[3]
    strand = ""
    if len(fields) > 4 and fields[4] in ["-", "+"]:
        strand = fields[4]

    return AnnotationRecord(chrom, start, end, title, strand)


def get_annotation_from_gtf_gff3_record(
    annotation_record_str, filename, line_number, file_type
):
    """
    Get the necessary info from a line extracted from a gtf or gff3 file,
    including strand and record title if present.

    * If the record is not for a gene, returns None
    to ignore it.
    * Prefers gene_name over gene_id for annotation title.

    Returns annotation record or None if not a valid gene record
    """
    fields = annotation_record_str.strip().split("\t")
    if len(fields) < 9:
        logger.error(
            "Annotation file {} contains invalid record at line {}".format(
                filename, line_number
            )
        )
        sys.exit(1)
    if fields[2] != "gene":
        return
    chrom = fields[0]
    start = int(fields[3])
    end = int(fields[4])

    # enforce start less than end
    if start > end:
        tmp = start
        start = end
        end = tmp
    title = ""
    gene_id = ""

    keyval_delim = "=" if file_type == ".gff3" else " "
    info_fields = fields[8].strip(" ;").split(";")
    for info_field in info_fields:
        if keyval_delim not in info_field:
            continue
        value_name, value = info_field.strip().split(keyval_delim)
        value = value.strip(' "').lower()
        value_name = value_name.strip().lower()

        if value_name == "gene_id":
            gene_id = value.strip('"')
        if value_name == "gene_name":
            title = value.strip('"')
    if title == "":
        title = gene_id

    strand = fields[6]
    if strand not in ["-", "+"]:
        logger.warning(
            "Invalid Strand field {} in annotation file {} at line {}".format(
                strand, filename, line_number
            )
        )
    return AnnotationRecord(chrom, start, end, title, strand)


def unpack_annotation_records(annotation_filenames):
    """
    Read annotation files (bed, gtf, or gff3) into a
    data structure that stores a dictionary with chromosome string
    as key, and an ordered dict of start coordinate to
    (end coordinate/strand) as value. Also optionally includes
    strand and annotation title. Returns a dict keyed by filename.
    """
    all_annotation_records = {}
    if type(annotation_filenames) is not list:
        annotation_filenames = [annotation_filenames]
    for annotation_filename in annotation_filenames:
        if annotation_filename is None or len(annotation_filename) == 0:
            return all_annotation_records
        annotations_fh = (
            gzip.open(annotation_filename, "rt", encoding="UTF-8")
            if is_gzipped(annotation_filename)
            else open(annotation_filename, "r")
        )

        # identify the correct function to get the relevant info from the file
        file_extension = os.path.splitext(os.path.basename(annotation_filename))[
            1
        ].lower()
        if file_extension == ".gz":
            file_extension = os.path.splitext(
                os.path.splitext(os.path.basename(annotation_filename))[0]
            )[1].lower()

        get_annotation = ""
        if file_extension == ".bed":
            get_annotation = get_annotation_from_bed_record
        elif file_extension == ".gff3" or file_extension == ".gtf":
            get_annotation = get_annotation_from_gtf_gff3_record
        else:
            logger.error(
                "Unrecognized annotation file type {} for file {}".format(
                    annotation_filename, file_extension
                )
            )
            sys.exit(1)

        annotation_records = {}
        for i, line in enumerate(annotations_fh):
            if len(line) == 0 or line[0] == "#":
                continue
            annotation = get_annotation(line, annotation_filename, i, file_extension)
            if annotation:
                if annotation.chrom not in annotation_records:
                    annotation_records[annotation.chrom] = []
                annotation_records[annotation.chrom].append(annotation)
        annotations_fh.close()

        # force sorted order
        for chrom in annotation_records:
            annotation_records[chrom].sort(key=lambda r: (r.end, r.start))
        annotation_name = os.path.basename(annotation_filename)
        all_annotation_records[annotation_name] = annotation_records
    return all_annotation_records
