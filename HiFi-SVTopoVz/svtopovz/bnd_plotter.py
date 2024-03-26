#!/usr/bin/env python

import re
import sys
import matplotlib.pyplot as plt
from svtopovz import bed_plotter
from svtopovz import chain_plotter
import numpy as np
import logging

plt.set_loglevel(level="warning")

logger = logging.getLogger(__name__)

ARROWHEAD_SIZE = 50
WINDOW_EDGE_SIZE = 20
MIN_CLUSTER_PRECISION = 10


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
    Converts a string with chr:pos-end to dictionary
    """
    try:
        window_components = region_window.split(":")
        chrom = window_components[0][0].upper() + window_components[0][1:]
        if "-" in window_components[1]:
            start, end = window_components[1].split("-")
        else:
            start, end = window_components[1], window_components[1]
        chrom = chrom.lower()
        start = int(start)
        end = int(end)
        return {"chrom": chrom, "start": start, "end": end}
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
            result.append(text.lower())
    return result


def get_full_region(region_info):
    """
    Gets the full genomic region or window from a list of region info entries
    """
    positions = []

    for entry in region_info:
        for coverage in entry["coverages"]:
            positions.append(coverage)
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
        print("ERROR: Region start is non-integer", file=sys.stderr)
        print(e, file=sys.stderr)
        sys.exit()
    try:
        end = int(end)
    except OverflowError as e:
        print("ERROR: Region end is non-integer", file=sys.stderr)
        print(e, file=sys.stderr)
        sys.exit()

    length = abs(end - start)
    len_str = ""
    if length <= 10_000:
        len_str = "{} bp".format(length)
    elif 10_000 < length < 1_000_000:
        # if greater than 10kb, round to nearest thousand and report as kb
        rounded_length = round(length / 1_000)
        len_str = "{} kbp".format(rounded_length)
    else:
        # if greater than 1mb, round to nearest million and report as mb
        rounded_length = round(length / 1_000_000)
        len_str = "{} mbp".format(rounded_length)
    return len_str


def plot_spanned_block(
    coverages: dict,
    plot_idx: float,
    window_size,
    orientation: bool,
    matched_order_count: int,
    ax,
):
    """
    Adds plot elements for one block of a complex SV with spanning alignments.
    If multiple blocks occur at the same sample order index, they must shrink to
    fit in the space provided. The plot index can therefore scale by this block's
    position in the sample order index, making it a float instead of integer.
    """
    end_coords = []
    start_coords = []
    sorted_coverages = sorted(list(coverages.keys()), key=natural_sort_key)
    if len(sorted_coverages) == 0:
        print(
            "Warning: no coverages found for plot index {}".format(plot_idx),
            file=sys.stderr,
        )
        return start_coords, end_coords

    begin_region_coordinates = get_coordinates(sorted_coverages[0])
    rectangle_start = begin_region_coordinates["start"]
    rectangle_end = get_coordinates(sorted_coverages[-1])["end"]

    # plot alignments coverage if this isn't the end
    # each entry is the end of that coverage depth so plot previous coverage each time
    prev_alignment_start = None
    prev_coverage = coverages[sorted_coverages[0]]
    cov_color = "black"
    for coverage_coord in sorted_coverages:
        curr_alignment_dict = get_coordinates(coverage_coord)
        if prev_alignment_start is not None:
            ax.hlines(
                y=plot_idx,
                xmin=prev_alignment_start,
                xmax=curr_alignment_dict["start"],
                linewidth=prev_coverage / matched_order_count,
                color=cov_color,
            )
            prev_coverage = coverages[coverage_coord]
        prev_alignment_start = curr_alignment_dict["start"]

    # plot rectangle around the block
    arrow_len = window_size / ARROWHEAD_SIZE
    rectangle_fill = False
    rectangle_color = "darkgrey"
    if matched_order_count > 1:
        rectangle_color = "pink"
    patch_halfwidth = 0.25 / matched_order_count  # based on the unit height of each row

    rectangle = plt.Polygon(
        [
            [rectangle_start, plot_idx + patch_halfwidth],
            [rectangle_start, plot_idx - patch_halfwidth],
            [rectangle_end, plot_idx - patch_halfwidth],
            [rectangle_end, plot_idx + patch_halfwidth],
            [rectangle_start, plot_idx + patch_halfwidth],
        ],
        fill=rectangle_fill,
        lw=0.5,
        color=rectangle_color,
    )
    ax.add_patch(rectangle)

    # default to forward annotation of ends
    end_coords = [rectangle_end, plot_idx]
    start_coords = [rectangle_start, plot_idx]

    # plot arrowhead for end of rectangle
    arrowhead = None
    if orientation == "+":
        arrowhead = plt.Polygon(
            [
                [rectangle_end, plot_idx + patch_halfwidth],
                [rectangle_end + arrow_len, plot_idx],
                [rectangle_end, plot_idx - patch_halfwidth],
                [rectangle_end, plot_idx + patch_halfwidth],
            ],
            color=rectangle_color,
            fill=rectangle_fill,
            lw=0.5,
        )
    elif orientation == "-":
        arrowhead = plt.Polygon(
            [
                [rectangle_start, plot_idx + patch_halfwidth],
                [rectangle_start - arrow_len, plot_idx],
                [rectangle_start, plot_idx - patch_halfwidth],
                [rectangle_start, plot_idx + patch_halfwidth],
            ],
            color=rectangle_color,
            fill=rectangle_fill,
            lw=0.5,
        )
        end_coords = [rectangle_start, plot_idx]
        start_coords = [rectangle_end, plot_idx]

    if arrowhead is not None:
        ax.add_patch(arrowhead)

    return start_coords, end_coords


def is_inverted(block):
    """
    checks if a block is inverted. Returns true if block start is greater than block start
    """
    return block["region"]["start"] > block["region"]["end"]


def generate_midpoints(pointa, pointb):
    """
    given two points, generates the 10 midpoints between them
    """
    midpoints = [pointa]
    for i in range(1, 10):
        x = pointa[0] + (pointb[0] - pointa[0]) * i / 10
        y = pointa[1] + (pointb[1] - pointa[1]) * i / 10
        midpoints.append((x, y))
    return midpoints


def get_unspanned_connections(current_block, spanned_block_ends):
    """
    From an unspanned connection, find the previous and next breakends it connects

    This is done by finding the closest end position to the start of the current block
    and the closest start to the end of the current block
    """
    potential_previous_connections = []
    potential_next_connections = []
    curr_end = current_block["region"]["end"]
    curr_start = current_block["region"]["start"]

    for end_key in spanned_block_ends["ends"]:
        spanned_genomic_coord = end_key[0]
        plot_idx = end_key[1]
        end_block = spanned_block_ends["ends"][end_key]

        if plot_idx - current_block["sample_order_index"] > 1:
            continue

        if abs(curr_start - spanned_genomic_coord) < MIN_CLUSTER_PRECISION:
            potential_previous_connections.append((end_block, plot_idx))
        if abs(curr_end - spanned_genomic_coord) < MIN_CLUSTER_PRECISION:
            potential_next_connections.append((end_block, plot_idx))

    for start_key in spanned_block_ends["starts"]:
        spanned_genomic_coord = start_key[0]
        plot_idx = start_key[1]
        start_block = spanned_block_ends["starts"][start_key]

        if plot_idx - current_block["sample_order_index"] > 1:
            continue

        if abs(curr_end - spanned_genomic_coord) < MIN_CLUSTER_PRECISION:
            potential_next_connections.append((start_block, plot_idx))
        if abs(curr_start - spanned_genomic_coord) < MIN_CLUSTER_PRECISION:
            potential_previous_connections.append((start_block, plot_idx))

    prev_connection = None
    next_connection = None
    prev_plot_idx = None
    next_plot_idx = None
    if len(potential_previous_connections) == 1:
        prev_connection, prev_plot_idx = potential_previous_connections[0]
    if len(potential_next_connections) == 1:
        next_connection, next_plot_idx = potential_next_connections[0]

    if next_connection and prev_connection is None:
        # if one of the connections was ambigous and the other was not, use
        # the unambigious to attempt to recover the ambiguous one
        # this can help if two blocks have the same start and different ends in a dup
        new_potential_prev_connections = []
        for potential_prev_connection, plot_idx in potential_previous_connections:
            if potential_prev_connection["region"] != next_connection["region"]:
                new_potential_prev_connections.append(
                    (potential_prev_connection, plot_idx)
                )
        if len(new_potential_prev_connections) == 1:
            prev_connection = new_potential_prev_connections[0]

    elif prev_connection and next_connection is None:
        # same as above but for next_connection recovery
        new_potential_next_connections = []
        for potential_next_connection, plot_idx in potential_next_connections:
            if potential_next_connection["region"] != prev_connection["region"]:
                new_potential_next_connections.append(
                    (potential_next_connection, plot_idx)
                )
        if len(new_potential_next_connections) == 1:
            next_connection, next_plot_idx = new_potential_next_connections[0]

    result = {
        "prev_block": prev_connection,
        "next_block": next_connection,
        "prev_plot_idx": prev_plot_idx,
        "next_plot_idx": next_plot_idx,
    }

    return result


def plot_unspanned_block(
    current_block,
    plot_idx: float,
    window_start: int,
    window_end: int,
    spanned_block_ends: dict,
    ax,
):
    """
    Adds plot elements for one block of a complex SV without spanning alignments
    """
    window_size = window_end - window_start
    arrow_len = window_size / ARROWHEAD_SIZE

    connection_info = get_unspanned_connections(current_block, spanned_block_ends)
    prev_break = connection_info["prev_block"]
    next_break = connection_info["prev_block"]
    prev_plot_idx = connection_info["prev_plot_idx"]
    next_plot_idx = connection_info["next_plot_idx"]

    arrow_start = current_block["region"]["start"]
    arrow_end = current_block["region"]["end"]
    logger.debug("Unspanned connection previous break: {}".format(prev_break))
    logger.debug("Unspanned connection next break: {}".format(next_break))
    if prev_break is None or next_break is None:
        return
    # margin for the arrowhead if uniquely present
    if prev_break is not None and type(prev_break) is not tuple:
        if prev_break["orientation"] == "+":
            arrow_start += arrow_len
        elif prev_break["orientation"] == "-":
            arrow_start -= arrow_len

    if prev_plot_idx is None:
        prev_plot_idx = plot_idx + 1
    if next_plot_idx is None:
        next_plot_idx = plot_idx - 1
    plot_idx = prev_plot_idx - (prev_plot_idx - next_plot_idx) / 2

    start = (arrow_start, prev_plot_idx)
    midpoint1 = (
        arrow_start,
        plot_idx,
    )

    midpoint2 = (
        arrow_end,
        plot_idx,
    )

    end = (arrow_end, next_plot_idx)
    arrow_coords = [
        start,
        midpoint1,
        midpoint2,
        end,
    ]

    x = []
    y = []
    for coord in arrow_coords:
        x.append(coord[0])
        y.append(coord[1])
    x = np.array(x)
    y = np.array(y)
    ax.plot(x, y, linestyle="dashed", color="gray", linewidth=0.75)


def plot_unspanned_blocks(
    region_info_by_sample_order,
    max_sample_idx,
    block_count,
    spanned_block_ends,
    coordinates,
    ax,
):
    """
    Go through all blocks in this regions and plot the ones that are not spanned.
    Use the stored spanned block ends to determine where the unspanned blocks start and end.
    """
    for sample_idx in range(max_sample_idx + 1):
        if sample_idx not in region_info_by_sample_order:
            continue
        for i, breakend in enumerate(region_info_by_sample_order[sample_idx]):
            matched_order_count = len(region_info_by_sample_order[sample_idx])
            coverages = breakend["coverages"]
            if len(coverages) > 0:
                # spanned, already done
                continue
            sample_idx = breakend["sample_order_index"]

            plot_idx = (
                block_count - sample_idx
            )  # plot_idx is from top of image to bottom
            if matched_order_count > 1:
                # subdivide plot index if multiple blocks in it
                plot_idx = (plot_idx - 1) + (i + 1) / matched_order_count

            plot_unspanned_block(
                breakend,
                plot_idx,
                coordinates["start"],
                coordinates["end"],
                spanned_block_ends,
                ax,
            )


def plot_spanned_blocks(
    region_info_by_sample_order, max_sample_idx, block_count, region_size, ax
):
    """
    Go through all blocks in this region and plot the ones that are spanned.
    Store their plot coordinates in a dict keyed by 'starts' or 'ends', then
    a dict of coordinate to block info.
    """
    spanned_block_ends = {
        "starts": {},
        "ends": {},
    }
    for sample_idx in range(max_sample_idx + 1):
        if sample_idx not in region_info_by_sample_order:
            continue
        for i, breakend in enumerate(region_info_by_sample_order[sample_idx]):
            matched_order_count = len(
                [
                    x
                    for x in region_info_by_sample_order[sample_idx]
                    if len(x["coverages"]) > 0
                ]
            )
            coverages = breakend["coverages"]
            if len(coverages) == 0:
                # unspanned, save for later
                continue
            sample_idx = breakend["sample_order_index"]
            orientation = breakend["orientation"]

            plot_idx = (
                block_count - sample_idx
            )  # plot_idx is from top of image to bottom
            if matched_order_count > 1:
                # subdivide plot index if multiple blocks in it
                plot_idx = (plot_idx - 1) + (i + 1) / matched_order_count

            # plots a block from the start of the alignment region to the end for a pair of
            # coordinates with spanning alignments, on the same level in the plot
            block_start_coords, block_end_coords = plot_spanned_block(
                coverages,
                plot_idx,
                region_size,
                orientation,
                matched_order_count,
                ax,
            )
            if len(block_start_coords) > 0:
                spanned_block_ends["starts"][tuple(block_start_coords)] = breakend
            if len(block_end_coords) > 0:
                spanned_block_ends["ends"][tuple(block_end_coords)] = breakend
    return spanned_block_ends


def plot_region(region_info, out_prefix, bed_records, image_type):
    """
    Generates a static image representation of one complex SV and any bed annotations
    """
    region_info_by_sample_order, max_sample_idx = order_by_sample_idx(region_info)
    coordinates = get_full_region(region_info)
    if coordinates is None:
        return
    region_size = coordinates["end"] - coordinates["start"]
    plot_height = len(region_info) + 1
    chain_height = min(3, int(plot_height / 2))

    fig, axes = plt.subplots(
        2, 1, gridspec_kw={"height_ratios": [plot_height, chain_height]}
    )
    fig.set_figheight(5)
    fig.set_figwidth(12)
    axes[0].set_xlim(
        coordinates["start"] - (region_size / WINDOW_EDGE_SIZE),
        coordinates["end"] + (region_size / WINDOW_EDGE_SIZE),
    )

    axes[0].set_ylim(-0.25, plot_height + 0.25)
    axes[1].set_ylim(0, 1)

    spanned_block_ends = plot_spanned_blocks(
        region_info_by_sample_order,
        max_sample_idx,
        len(region_info),
        region_size,
        axes[0],
    )
    plot_unspanned_blocks(
        region_info_by_sample_order,
        max_sample_idx,
        len(region_info),
        spanned_block_ends,
        coordinates,
        axes[0],
    )

    chain_plotted = chain_plotter.plot_chain_representation(
        region_info_by_sample_order,
        coordinates,
        plot_height,
        axes,
    )
    if chain_plotted:
        # remove ticks for chain plot
        axes[1].set_xticks([])
        axes[1].set_yticks([])
        for spine in axes[1].spines.values():
            spine.set_visible(False)
        plt.subplots_adjust(hspace=0.75)
    else:
        fig.delaxes(axes[1])

    # configure ticks for bnd plot
    axes[0].ticklabel_format(style="plain")
    axes[0].set_yticks([])
    xtick_locations = axes[0].get_xticks()
    xtick_labels = ["{:,}".format(int(l)) for l in xtick_locations]
    axes[0].set_xticks(xtick_locations, xtick_labels, rotation=30)

    plt.suptitle(
        "{}:{}-{}".format(
            coordinates["chrom"], coordinates["start"], coordinates["end"]
        )
    )
    axes[0].set_title(
        round_genomic_coords_to_str(coordinates["start"], coordinates["end"]),
        fontsize=10,
    )
    bed_plotter.annotate(bed_records, coordinates, axes[0])

    figname = "{}_{}-{}-{}.{}".format(
        out_prefix,
        coordinates["chrom"],
        coordinates["start"],
        coordinates["end"],
        image_type,
    )
    plt.savefig(figname, dpi=400)
    plt.close()
