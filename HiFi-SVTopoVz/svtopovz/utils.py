#!/usr/bin/env python

import re
import sys
import numpy as np
import logging
from collections import namedtuple

logger = logging.getLogger(__name__)

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
Region = namedtuple("Region", ["chrom", "start", "end", "size"])


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
        chrom = window_components[0][0].upper() + window_components[0][1:]
        if "-" in window_components[1]:
            start, end = window_components[1].split("-")
        else:
            start, end = window_components[1], window_components[1]
        chrom = chrom.upper()
        start = int(start)
        end = int(end)
        if start < end:
            start, end = end, start
        return Region(chrom, start, end, end - start)
    except ValueError:
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
        start["end"] = end["end"]
        return start


def round_genomic_coords_to_str(start, end):
    try:
        start = int(start)
    except OverflowError as e:
        logger.error("ERROR: Region start is non-integer")
        logger.error(e)
        sys.exit()
    try:
        end = int(end)
    except OverflowError as e:
        logger.error("ERROR: Region end is non-integer")
        logger.error(e)
        sys.exit()

    length = abs(end - start)
    len_str = ""
    if length < 1_000:
        len_str = "{} bp".format(length)
    elif 1_000 <= length <= 10_000:
        # if greater than 1kb, round to nearest hundred and report as kb
        rounded_length = int(round(length / 1_00) / 10)
        len_str = "{} kbp".format(rounded_length)
    elif 10_000 < length < 1_000_000:
        # if greater than 10kb, round to nearest thousand and report as kb
        rounded_length = round(length / 1_000)
        len_str = "{} kbp".format(rounded_length)
    else:
        # if greater than 1mb, round to nearest million and report as mb
        rounded_length = round(length / 1_000_000)
        len_str = "{} mbp".format(rounded_length)
    return len_str
