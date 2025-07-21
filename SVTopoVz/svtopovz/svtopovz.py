#!/usr/bin/env python
from collections import namedtuple

import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import logging
from glob import glob
from svtopovz.utils import *
from svtopovz.bnd_plotter import plot_breaks
from svtopovz.chain_plotter import plot_chain_representation
from svtopovz.annotation_plotter import annotate
from svtopovz.page_builder import build_review_page
from .__init__ import __version__
from tqdm import tqdm
from collections import namedtuple

logger = logging.getLogger(__name__)
logging.getLogger("matplotlib.font_manager").disabled = True
ARROWHEAD_SIZE = 50
MakePlotArgs = namedtuple(
    "MakePlotArgs",
    [
        "event_info",
        "gene_annotation_records",
        "bed_records",
        "max_gap_size_mb",
        "prefix",
        "svtopo_dir",
        "image_type",
    ],
)


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
    positions = get_positions_in_window(region_info)
    windows_by_chrom = get_windows_by_chromosome(positions, max_size)
    windows = get_windows(windows_by_chrom, positions)
    return windows


def get_positions_in_window(region_info):
    """
    Given region info, extract all break positions from the event
    """
    positions = {}
    for entry in region_info:
        start_chrom = entry["region"]["start_chrom"]
        end_chrom = entry["region"]["end_chrom"]
        start = entry["region"]["start"]
        end = entry["region"]["end"]

        # Initialize chromosome lists if needed
        if start_chrom not in positions:
            positions[start_chrom] = set()
        if end_chrom not in positions:
            positions[end_chrom] = set()

        # Add positions to sets to avoid duplicates
        positions[start_chrom].add(start)
        positions[end_chrom].add(end)

    # Convert sets to sorted lists
    for chrom in positions:
        positions[chrom] = sorted(list(positions[chrom]))

    return positions


def get_windows_by_chromosome(positions, max_size):
    """
    Given positions in the event, return dict of windows of at most max_size that
    include all of the positions keyed by chrom
    """
    windows_by_chrom = {}
    for chrom in positions:
        if chrom not in windows_by_chrom and len(positions[chrom]) > 0:
            windows_by_chrom[chrom] = []
        chrom_positions = sorted(list(set(positions[chrom])))

        if len(chrom_positions) == 0:
            continue

        if len(chrom_positions) == 1:
            # Single position - create a window around it
            windows_by_chrom[chrom].append(
                Region(chrom, chrom_positions[0], chrom_positions[0], 0)
            )
            continue

        window_start_idx = 0
        for i, end_pos in enumerate(chrom_positions):
            start_pos = chrom_positions[window_start_idx]
            window_size = end_pos - start_pos
            if window_size >= max_size:
                prev_end_pos = (
                    chrom_positions[i - 1]
                    if i > window_start_idx
                    else start_pos + SPLIT_REGION_PADDING
                )
                window_size = prev_end_pos - start_pos
                windows_by_chrom[chrom].append(
                    Region(chrom, start_pos, prev_end_pos, window_size)
                )
                window_start_idx = i
        start_pos = chrom_positions[window_start_idx]
        end_pos = chrom_positions[-1]
        window_size = end_pos - start_pos
        windows_by_chrom[chrom].append(Region(chrom, start_pos, end_pos, window_size))
    return windows_by_chrom


def get_windows(windows_by_chrom, positions):
    """
    Given a dict of windows by their chromosome and all positions in an event
    return all windows split at chromosome and at the maximum window size.
    """
    windows = []
    for chrom in windows_by_chrom:
        for entry in windows_by_chrom[chrom]:
            windows.append(entry)
    windows = sorted(windows, key=lambda w: [w.chrom, w.start, w.end])
    if len(windows) > 0:
        # Add padding to first window
        first_window = windows[0]
        windows[0] = Region(
            first_window.chrom,
            first_window.start - SPLIT_REGION_PADDING,
            first_window.end,
            first_window.end - (first_window.start - SPLIT_REGION_PADDING),
        )

        # Add padding to last window
        last_window = windows[-1]
        windows[-1] = Region(
            last_window.chrom,
            last_window.start,
            last_window.end + SPLIT_REGION_PADDING,
            (last_window.end + SPLIT_REGION_PADDING) - last_window.start,
        )
    elif len(windows) == 0 and positions:
        # Get the first chromosome's positions and convert to integers
        chrom = list(positions.keys())[0]  # Get the first chromosome
        chrom_positions = [int(pos) for pos in positions[chrom]]
        start = min(chrom_positions)
        end = max(chrom_positions)
        # Add padding to both sides since it's a single window
        windows.append(
            Region(
                chrom,
                start - SPLIT_REGION_PADDING,
                end + SPLIT_REGION_PADDING,
                (end + SPLIT_REGION_PADDING) - (start - SPLIT_REGION_PADDING),
            )
        )
    return windows


def create_complex_sv_image(event_info, gene_records, bed_records, max_window_size):
    """
    Create an image of a complex SV, using the JSON record for the event.
    This function creates the pyplot.Figure object with the correct subplots,
    including the definition and control of split view, and calls functions
    from the bnd_plotter, chain_plotter, and bed_plotter to fill in the
    image.
    """
    if not event_info:
        return []
    region_info_by_sample_order, max_sample_idx = order_by_sample_idx(event_info)
    windows = choose_window_coordinates(event_info, max_window_size)
    if windows is None:
        return []
    annotation_region_height = len(bed_records) / 1.5
    plot_height = len(event_info) + 2 + annotation_region_height
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
            annotation_region_height,
            ax,
        )
        annotation_level = 0
        for gene_filename in gene_records:
            annotate(
                gene_filename,
                gene_records[gene_filename],
                window,
                axis_dict[BREAKS][-1],
                annotation_level,
                i,
            )
            annotation_level += 1
        for bed_record_group in bed_records:
            annotate(
                bed_record_group,
                bed_records[bed_record_group],
                window,
                axis_dict[BREAKS][-1],
                annotation_level,
                i,
            )
            annotation_level += 1
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
    plt.subplots_adjust(hspace=0.35, right=0.9)

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


def is_filtered_out(count, event_info, include_simple_breakpoints):
    """
    Compare the event info against filtering options and
    definitions of filterable events to determine if
    this event should be plotted.
    """
    # Handle empty events first
    if len(event_info) == 0:
        return True

    orientations = [x["orientation"] for x in event_info]
    entry_types = [len(x["coverages"]) > 0 for x in event_info]
    chroms = set([x["region"]["start_chrom"] for x in event_info]).union(
        set([x["region"]["end_chrom"] for x in event_info])
    )
    is_simple_del = (
        orientations == ["+", "+", "+"]
        and entry_types == [True, False, True]
        and len(chroms) == 1
    )

    is_simple_dup = (
        orientations == ["+", "-", "+"]
        and entry_types == [True, False, True]
        and len(chroms) == 1
    )
    is_nonreciprocal_translocation = (
        entry_types == [True, False, True] and len(chroms) == 2
    )
    is_simple_breakpoint = entry_types == [True, False, True]

    event_type = "unknown"
    if is_simple_del:
        event_type = "Deletion"
    elif is_simple_dup:
        event_type = "Duplication"
    elif is_nonreciprocal_translocation:
        event_type = "Nonreciprocal translocation"
    elif len(event_info) == 1:
        event_type = "Single-ended BND"
        is_simple_breakpoint = True

    if not include_simple_breakpoints and is_simple_breakpoint:
        logger.debug("Skipped event #{}. Event type: {}".format(count, event_type))
        return True
    return False


# get the region-based name for the image, omitting the path and extension
def get_image_name(event_info):
    image_regions = {}
    for block in event_info:
        if len(block["coverages"]) != 0:
            chrom = block["region"]["start_chrom"]
            if chrom not in image_regions:
                image_regions[chrom] = []
            image_regions[chrom].append(block["region"]["start"])
            image_regions[chrom].append(block["region"]["end"])

    chroms = sorted(image_regions.keys())
    simplified_regions = []
    for chrom in chroms:
        start = min(image_regions[chrom])
        end = max(image_regions[chrom])
        simplified_regions.append("{}-{}-{}".format(chrom, start, end))
    return "__".join(simplified_regions)


def validate_sv_info(sv_info):
    """
    Make sure that the expected fields exist in the input json object
    """
    assert "svtopo_version" in sv_info

    input_version = sv_info["svtopo_version"]
    if __version__ != input_version:
        logger.warning(
            "Input data version {} does not match this version {}".format(
                input_version, __version__
            )
        )

    assert "event_graphs" in sv_info


def make_plot(plot_args: MakePlotArgs):
    """
    Create and save a single image using a MakePlotArgs tuple as input
    """
    figname_splits = create_complex_sv_image(
        plot_args.event_info,
        plot_args.gene_annotation_records,
        plot_args.bed_records,
        plot_args.max_gap_size_mb,
    )
    plt.suptitle(plot_args.prefix)
    if figname_splits:
        figname = os.path.join(
            plot_args.svtopo_dir,
            "images",
            "{}_{}.{}".format(
                plot_args.prefix,
                get_image_name(plot_args.event_info),
                plot_args.image_type,
            ),
        )

        plt.subplots_adjust(
            left=0.01,  # Reduce left margin (default is ~0.125)
            bottom=0.0,  # Reduce bottom margin (default is ~0.11)
            right=0.9,  # Keep right margin for legend space
            top=0.9,  # Keep top margin for title space
        )

        plt.savefig(figname, dpi=400)
        logger.debug("Saved image {}".format(figname))
        plt.close()


def make_plots(
    jsons,
    gene_annotation_records,
    bed_records,
    include_simple_breakpoints,
    max_gap_size_mb,
    svtopo_dir,
    image_type,
):
    """
    Generates all plots, with their annotations
    """
    os.makedirs(os.path.join(svtopo_dir, "images"), exist_ok=True)
    make_plot_args = []
    for json in jsons:
        unpacked = unpack_json(json)
        if unpacked is None:
            return
        sv_info, prefix = unpacked
        validate_sv_info(sv_info)
        for count, event_info in enumerate(sv_info["event_graphs"]):
            if not is_filtered_out(count, event_info, include_simple_breakpoints):
                make_plot_args.append(
                    MakePlotArgs(
                        event_info,
                        gene_annotation_records,
                        bed_records,
                        max_gap_size_mb * MEGABASE,
                        prefix,
                        svtopo_dir,
                        image_type,
                    )
                )
    for plot_args in tqdm(make_plot_args, desc="Creating images:"):
        make_plot(plot_args)


def svtopovz(args):
    jsons = glob(os.path.join(args.svtopo_dir, "*json"))
    jsons.extend(glob(os.path.join(args.svtopo_dir, "*json.gz")))
    gene_annotation_records = unpack_annotation_records(args.genes)
    bed_records = unpack_annotation_records(args.annotation_bed)

    make_plots(
        jsons,
        gene_annotation_records,
        bed_records,
        args.include_simple_breakpoints,
        args.max_gap_size_mb,
        args.svtopo_dir,
        args.image_type,
    )
    build_review_page(args)
