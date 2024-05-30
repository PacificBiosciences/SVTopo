#!/usr/bin/env python
from collections import namedtuple

import gzip
import json
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import logging
from svtopovz.utils import *
from svtopovz.bnd_plotter import plot_breaks
from svtopovz.chain_plotter import plot_chain_representation
from svtopovz.bed_plotter import annotate_bed_records

logger = logging.getLogger(__name__)
ARROWHEAD_SIZE = 50


def is_gzipped(putative_zipfile):
    """
    Check if file is zipped. Assumes file exists
    """
    with open(putative_zipfile, "rb") as filehandle:
        id_bytes = filehandle.read(2)
        return id_bytes == b"\x1f\x8b"


def unpack_json(json_filename):
    if not (os.path.exists(json_filename) and os.path.isfile(json_filename)):
        logger.error(" {} does not exist".format(json_filename))
        sys.exit()
    if json_filename.endswith(".gz"):
        if not is_gzipped(json_filename):
            logger.error("{} is identified as gzipped but is not".format(json_filename))
            sys.exit()
        with gzip.open(json_filename, "rt", encoding="UTF-8") as json_fh:
            try:
                return json.load(json_fh)
            except json.decoder.JSONDecodeError:
                logger.error(" {} is empty or misformatted".format(json_filename))
                sys.exit()
    else:
        with open(json_filename, "r") as json_fh:
            try:
                return json.load(json_fh)
            except json.decoder.JSONDecodeError:
                logger.error(" {} is empty or misformatted".format(json_filename))
                sys.exit()


def unpack_bed_records(bed_filename):
    """
    Read a bed file into a data structure that stores a dictionary with chromosome string as
    key, and an ordered dict of start coordinate to (end coordinate/strand) as value
    """
    if bed_filename is None or len(bed_filename) == 0:
        return
    bed_fh = (
        gzip.open(bed_filename, "rt", encoding="UTF-8")
        if bed_filename.endswith(".gz")
        else open(bed_filename, "r")
    )
    bed_record = namedtuple(
        "Record",
        [
            "start",
            "end",
            "title",
            "strand",
        ],
    )
    bed_records = {}
    for i, line in enumerate(bed_fh):
        if len(line) == 0 or line[0] == "#":
            continue
        fields = line.strip().split()
        if len(fields) < 3:
            logger.error(
                "BED file {} contains invalid record at line {}".format(bed_filename, i)
            )
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
        if chrom not in bed_records:
            bed_records[chrom] = []
        # sort the entries by the end position, than by the start position
        bed_records[chrom].append(bed_record(start, end, title, strand))
    bed_fh.close()

    # force sorted order
    for chrom in bed_records:
        bed_records[chrom].sort(key=lambda r: (r.end, r.start))
    return bed_records


def choose_window_coordinates(region_info, max_size):
    """
    Given a list of block coordinate dicts from JSON, select the
    coordinates that will be used for each split of coordinates,
    such that the max unspanned gap size is maintained. Return an ordered
    list of Regions for genomic windows. A different chromosome is
    always assumed to be greater than max_size away.
    """
    if max_size == 0:
        max_size = BIG_NUMBER
    positions = {}
    for entry in region_info:
        start_chrom = entry["region"]["start_chrom"]
        end_chrom = entry["region"]["end_chrom"]
        if start_chrom not in positions:
            positions[start_chrom] = []
        if end_chrom not in positions:
            positions[end_chrom] = []

        start = entry["region"]["start"]
        end = entry["region"]["end"]
        if start_chrom == end_chrom:
            positions[start_chrom].append(start)
            positions[start_chrom].append(end)
        else:
            positions[start_chrom].append(start)
            positions[end_chrom].append(end)

    windows_by_chrom = {}
    for chrom in positions:
        if chrom not in windows_by_chrom:
            windows_by_chrom[chrom] = []
        chrom_positions = sorted(positions[chrom])
        prev_start = 0
        prev_end = 0
        for pos in chrom_positions:
            if prev_start == 0:
                prev_start = pos
            elif prev_end > prev_start and (pos - prev_end) > max_size:
                current_dist = pos - prev_end
                padding = min(int(current_dist / 4), SPLIT_REGION_PADDING)
                windows_by_chrom[chrom].append(
                    Region(
                        chrom,
                        prev_start,
                        prev_end + padding,
                        (prev_end - prev_start) + padding,
                    )
                )
                prev_start = pos - padding
            prev_end = pos
        if len(windows_by_chrom[chrom]) > 0:
            first_region = windows_by_chrom[chrom][-1]
            windows_by_chrom[chrom][-1] = Region(
                chrom,
                first_region.start,
                prev_end,
                prev_end - first_region.start,
            )
        else:
            windows_by_chrom[chrom].append(
                Region(chrom, prev_start, prev_end, prev_end - prev_start)
            )
    windows = []
    for chrom in windows_by_chrom:
        for entry in windows_by_chrom[chrom]:
            windows.append(entry)
    windows = sorted(windows, key=lambda w: [w.chrom, w.start, w.end])
    if len(windows) > 0:
        start_window = windows[0]
        windows_by_chrom[0] = Region(
            start_window.chrom,
            start_window.start - SPLIT_REGION_PADDING,
            start_window.end,
            start_window.size + SPLIT_REGION_PADDING,
        )
        end_window = windows[-1]
        windows[-1] = Region(
            end_window.chrom,
            end_window.start,
            end_window.end + SPLIT_REGION_PADDING,
            end_window.size + SPLIT_REGION_PADDING,
        )
    elif len(windows) == 0:
        start = min(positions)
        end = max(positions)
        windows.append(Region(chrom, start, end, end - start))

    return windows


def create_complex_sv_image(event_info, bed_records, max_window_size):
    """
    Create an image of a complex SV, using the JSON record for the event.
    This function creates the pyplot.Figure object with the correct subplots,
    including the definition and control of split view, and calls functions
    from the bnd_plotter, chain_plotter, and bed_plotter to fill in the
    image.
    """
    region_info_by_sample_order, max_sample_idx = order_by_sample_idx(event_info)
    windows = choose_window_coordinates(event_info, max_window_size)
    if windows is None:
        return []
    plot_height = len(event_info) + 1
    chain_height = min(3, int(plot_height / 2))

    fig = plt.figure(figsize=(IMAGE_WIDTH, IMAGE_HEIGHT))
    height_ratios = [plot_height, chain_height]
    grid_specification = gridspec.GridSpec(
        ncols=len(windows), nrows=2, figure=fig, height_ratios=height_ratios
    )
    axis_dict = {
        BREAKS: [],
        CHAINS: None,
    }
    ##########################################################
    # set up break plot panes with x/y limits and titles,
    # then add break plots
    figname_splits = []
    for i, window in enumerate(windows):
        figname_splits.append("{0.chrom}_{0.start}_{0.end}".format(window))
        ax = fig.add_subplot(grid_specification[0, i])
        axis_dict[BREAKS].append(ax)

        ax.set_ylim(-0.25, plot_height + 0.5)
        plot_breaks(
            region_info_by_sample_order,
            max_sample_idx,
            window,
            len(event_info),
            len(windows),
            ax,
        )
        annotate_bed_records(bed_records, window, axis_dict[BREAKS][-1])
        box = [x for x in ax.spines.values()]
        if i > 0:
            box[0].set_linewidth(0.1)  # interior right spline
        if i < (len(windows) - 1):
            box[1].set_linewidth(0.1)  # interior left spline

    ##########################################################
    # add chain plots
    axis_dict[CHAINS] = fig.add_subplot(grid_specification[1, :])
    axis_dict[CHAINS].set_ylim(0, 1)
    axis_dict[CHAINS].set_xticks([])
    axis_dict[CHAINS].set_yticks([])
    for spine in axis_dict[CHAINS].spines.values():
        spine.set_visible(False)
    plt.subplots_adjust(hspace=0.75, right=0.9)

    chain_plotted = plot_chain_representation(
        region_info_by_sample_order,
        windows,
        axis_dict,
    )

    if not chain_plotted:
        logger.debug(
            "Failed to construct chain plot for {}:{}-{}".format(
                window.chrom, window.start, window.end
            )
        )
        plt.close()
        return
    return figname_splits


def svtopovz(args):
    sv_info = unpack_json(args.json)
    bed_records = unpack_bed_records(args.bed)
    for count, event_info in enumerate(sv_info):
        orientations = [x["orientation"] for x in event_info]
        entry_types = [len(x["coverages"]) > 0 for x in event_info]
        if len(event_info) == 0:
            logger.debug("Skipped event #{}: empty entry".format(count))
        elif len(event_info) == 1:
            logger.debug("Skipped event #{}: single-ended BND".format(count))
        elif (
            args.ignore_simple_dels
            and orientations == ["+", "+", "+"]
            and entry_types == [True, False, True]
        ):
            logger.debug("Skipped event #{}: simple DEL".format(count))
        elif (
            args.ignore_simple_dups
            and orientations == ["+", "-", "+"]
            and entry_types == [True, False, True]
        ):
            logger.debug("Skipped event #{}: simple DUP".format(count))
        else:
            figname_splits = create_complex_sv_image(
                event_info,
                bed_records,
                args.max_gap_size_mb * MEGABASE,
            )
            figname = "{}-{}.{}".format(
                args.out_prefix,
                "-".join(figname_splits),
                args.image_type,
            )

            plt.savefig(figname, dpi=400)
            plt.close()
            logger.debug("Finished event #{}".format(count))
