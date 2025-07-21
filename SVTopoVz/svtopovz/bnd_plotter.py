#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
from svtopovz.utils import *
import numpy as np
from copy import deepcopy
import logging
from svtopovz.utils import MAX_ALLOWED_COVERAGE

logger = logging.getLogger(__name__)


def plot_spanned_block(
    coverages: dict,
    plot_idx: float,
    window,
    orientation: bool,
    matched_order_count: int,
    sample_idx,
    annotation_region_height,
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
    # if the coverage is high, scale it down for better plotting
    scaled_coverages = scale_coverages(coverages)
    sorted_coverages = sorted(list(scaled_coverages.keys()), key=natural_sort_key)
    if len(sorted_coverages) == 0:
        logger.warning(
            "Warning: no coverages found for plot index {}".format(plot_idx),
        )
        return start_coords, end_coords

    begin_region_coordinates = get_coordinates(sorted_coverages[0])
    rectangle_start = begin_region_coordinates.start
    rectangle_end = get_coordinates(sorted_coverages[-1]).end

    if begin_region_coordinates.chrom.upper() != window.chrom.upper():
        return start_coords, end_coords
    starts_in_window = window.start <= rectangle_start <= window.end
    ends_in_window = window.start <= rectangle_end <= window.end
    if not starts_in_window and ends_in_window:
        return start_coords, end_coords

    plot_coverages(
        ax,
        scaled_coverages,
        sorted_coverages,
        matched_order_count,
        window.chrom,
        plot_idx + annotation_region_height,
    )
    plot_rectangle(
        ax,
        matched_order_count,
        rectangle_start,
        rectangle_end,
        plot_idx + annotation_region_height,
    )

    start_coords, end_coords = plot_arrowhead(
        ax,
        window.size,
        matched_order_count,
        rectangle_start,
        rectangle_end,
        orientation,
        sample_idx,
        plot_idx + annotation_region_height,
    )

    return start_coords, end_coords


def plot_coverages(
    ax,
    scaled_coverages,
    sorted_coverages,
    matched_order_count,
    window_chrom,
    y_value,
):
    """
    Plot alignments coverage if this isn't the end.
    Each entry is the end of that coverage depth so plot previous coverage each time
    """
    prev_alignment_start = None
    prev_coverage = scaled_coverages[sorted_coverages[0]]
    for coverage_coord in sorted_coverages:
        curr_alignment = get_coordinates(coverage_coord)
        if curr_alignment.chrom.upper() != window_chrom.upper():
            continue
        if prev_alignment_start is not None:
            ax.hlines(
                y=y_value,
                xmin=prev_alignment_start,
                xmax=curr_alignment.start,
                linewidth=prev_coverage / matched_order_count,
                color=COVERAGE_COLOR,
            )
            prev_coverage = scaled_coverages[coverage_coord]
        prev_alignment_start = curr_alignment.start


def plot_rectangle(
    ax,
    matched_order_count,
    rectangle_start,
    rectangle_end,
    y_value,
):
    """
    Plot the rectangle around a spanned block
    """
    rectangle_fill = False
    rectangle_color = "darkgrey"
    if matched_order_count > 1:
        rectangle_color = "pink"
    patch_halfwidth = 0.25 / matched_order_count  # based on the unit height of each row

    rectangle = plt.Polygon(
        [
            [rectangle_start, y_value + patch_halfwidth],
            [rectangle_start, y_value - patch_halfwidth],
            [rectangle_end, y_value - patch_halfwidth],
            [rectangle_end, y_value + patch_halfwidth],
            [rectangle_start, y_value + patch_halfwidth],
        ],
        fill=rectangle_fill,
        lw=0.5,
        color=rectangle_color,
    )
    ax.add_patch(rectangle)


def plot_arrowhead(
    ax,
    window_size,
    matched_order_count,
    rectangle_start,
    rectangle_end,
    orientation,
    sample_idx,
    y_value,
):
    """
    Plot the arrowhead that shows directionality of the block
    """
    # default to forward annotation of ends
    end_coords = [rectangle_end, sample_idx]
    start_coords = [rectangle_start, sample_idx]
    patch_halfwidth = 0.25 / matched_order_count  # based on the unit height of each row
    arrow_len = window_size / ARROWHEAD_SIZE
    rectangle_color = UNIQUE_BOX_COLOR if matched_order_count == 1 else AMBIG_BOX_COLOR

    # plot arrowhead for end of rectangle
    arrowhead = None
    if orientation == "+":
        arrowhead = plt.Polygon(
            [
                [rectangle_end, y_value + patch_halfwidth],
                [rectangle_end + arrow_len, y_value],
                [rectangle_end, y_value - patch_halfwidth],
                [rectangle_end, y_value + patch_halfwidth],
            ],
            color=rectangle_color,
            fill=False,
            lw=0.5,
        )
    elif orientation == "-":
        arrowhead = plt.Polygon(
            [
                [rectangle_start, y_value + patch_halfwidth],
                [rectangle_start - arrow_len, y_value],
                [rectangle_start, y_value - patch_halfwidth],
                [rectangle_start, y_value + patch_halfwidth],
            ],
            color=rectangle_color,
            fill=False,
            lw=0.5,
        )
        start_coords = [rectangle_end, sample_idx]
        end_coords = [rectangle_start, sample_idx]

    if arrowhead is not None:
        ax.add_patch(arrowhead)

    return start_coords, end_coords


def scale_coverages(coverages):
    """
    Scale high-coverages down if above the maximum
    """
    max_coverage = max(coverages.values())
    coverage_scaling_factor = max_coverage / MAX_ALLOWED_COVERAGE
    if max_coverage > MAX_ALLOWED_COVERAGE:
        for coord in coverages:
            coverages[coord] /= coverage_scaling_factor
    return coverages


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
    (
        potential_previous_connections,
        potential_next_connections,
    ) = get_unspanned_connection_candidates(current_block, spanned_block_ends)
    prev_connection, prev_sample_idx = refine_unspanned_connections(
        potential_previous_connections, current_block["sample_order_index"]
    )
    next_connection, next_sample_idx = refine_unspanned_connections(
        potential_next_connections, current_block["sample_order_index"]
    )

    if next_connection and prev_connection is None:
        # if one of the connections was ambigous and the other was not, use
        # the unambigious to attempt to recover the ambiguous one
        # this can help if two blocks have the same start and different ends in a dup
        new_potential_prev_connections = []
        for potential_prev_connection, sample_idx in potential_previous_connections:
            if potential_prev_connection["region"] != next_connection["region"]:
                new_potential_prev_connections.append(
                    (potential_prev_connection, sample_idx)
                )
        if len(new_potential_prev_connections) == 1:
            prev_connection = new_potential_prev_connections[0]

    elif prev_connection and next_connection is None:
        # same as above but for next_connection recovery
        new_potential_next_connections = []
        for potential_next_connection, sample_idx in potential_next_connections:
            if potential_next_connection["region"] != prev_connection["region"]:
                new_potential_next_connections.append(
                    (potential_next_connection, sample_idx)
                )
        if len(new_potential_next_connections) == 1:
            next_connection, next_sample_idx = new_potential_next_connections[0]

    result = {
        "prev_block": prev_connection,
        "next_block": next_connection,
        "prev_sample_idx": prev_sample_idx,
        "next_sample_idx": next_sample_idx,
    }
    return result


def get_unspanned_connection_candidates(current_block, spanned_block_ends):
    """
    Find spanned block starts and ends near the target block start and end
    """
    potential_previous_connections = []
    potential_next_connections = []
    curr_end = current_block["region"]["end"]
    curr_start = current_block["region"]["start"]
    if current_block["orientation"] == "-":
        curr_end, curr_start = curr_start, curr_end

    for end_key in spanned_block_ends["ends"]:
        spanned_genomic_coord = end_key[0]
        sample_idx = end_key[1]
        end_block = spanned_block_ends["ends"][end_key]

        if (sample_idx - current_block["sample_order_index"]) > 2:
            continue

        if abs(curr_start - spanned_genomic_coord) < MIN_CLUSTER_PRECISION:
            potential_previous_connections.append((end_block, sample_idx))
        if abs(curr_end - spanned_genomic_coord) < MIN_CLUSTER_PRECISION:
            potential_next_connections.append((end_block, sample_idx))

    for start_key in spanned_block_ends["starts"]:
        spanned_genomic_coord = start_key[0]
        sample_idx = start_key[1]
        start_block = spanned_block_ends["starts"][start_key]

        if (sample_idx - current_block["sample_order_index"]) > 2:
            continue

        if abs(curr_end - spanned_genomic_coord) < MIN_CLUSTER_PRECISION:
            potential_next_connections.append((start_block, sample_idx))
        if abs(curr_start - spanned_genomic_coord) < MIN_CLUSTER_PRECISION:
            potential_previous_connections.append((start_block, sample_idx))

    return potential_previous_connections, potential_next_connections


def refine_unspanned_connections(potential_connections, current_sample_index):
    """
    Given a set of connections that could match the start or end of
    the current block, refine down to 0-1 that are closest in sample space.
    """
    selected_connection = None
    selected_sample_idx = None

    min_order_dist = None
    for (
        potential_prev_connection,
        potential_prev_sample_idx,
    ) in potential_connections:
        current_dist = abs(current_sample_index - potential_prev_sample_idx)
        if min_order_dist is None or current_dist < min_order_dist:
            min_order_dist = current_dist
            selected_connection, selected_sample_idx = (
                potential_prev_connection,
                potential_prev_sample_idx,
            )
    return selected_connection, selected_sample_idx


def plot_unspanned_block(
    current_block,
    plot_idx: float,
    window,
    spanned_block_ends: dict,
    annotation_region_height,
    ax,
):
    """
    Adds plot elements for one block of a complex SV without spanning alignments
    """
    arrow_len = window.size / ARROWHEAD_SIZE

    connection_info = get_unspanned_connections(current_block, spanned_block_ends)
    prev_break = connection_info["prev_block"]
    next_break = connection_info["next_block"]
    prev_sample_idx = connection_info["prev_sample_idx"]
    next_sample_idx = connection_info["next_sample_idx"]

    current_block = reset_coordinates(window, deepcopy(current_block))
    arrow_start = current_block["region"]["start"]
    arrow_end = current_block["region"]["end"]

    if current_block["orientation"] == "-":
        arrow_end, arrow_start = arrow_start, arrow_end
    logger.debug("Unspanned connection from: {}".format(current_block))
    logger.debug("Previous break: {}".format(prev_break))
    logger.debug("Next break: {}".format(next_break))

    # margin for the arrowhead if uniquely present
    if (
        prev_break is not None
        and type(prev_break) is not tuple
        and current_block["orientation"] != "-"
    ):
        if prev_break["orientation"] == "+":
            arrow_start += arrow_len
        elif prev_break["orientation"] == "-":
            arrow_start -= arrow_len
    if prev_sample_idx is None:
        prev_sample_idx = plot_idx + 1
    if next_sample_idx is None:
        next_sample_idx = plot_idx - 1

    # the only options for the y-axis positions are matched with the current block or before and after
    prev_plot_index = plot_idx + 1
    next_plot_index = plot_idx - 1
    if prev_sample_idx == current_block["sample_order_index"]:
        plot_idx -= 0.5
        prev_plot_index -= 0.5
        next_plot_index += 0.5
    elif next_sample_idx == current_block["sample_order_index"]:
        next_plot_index += 1
        plot_idx += 0.5

    start = (arrow_start, prev_plot_index)
    midpoint1 = (
        arrow_start,
        plot_idx,
    )

    midpoint2 = (
        arrow_end,
        plot_idx,
    )

    end = (arrow_end, next_plot_index)
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
        y.append(coord[1] + annotation_region_height)
    x = np.array(x)
    y = np.array(y)
    ax.plot(x, y, linestyle="dashed", color="gray", linewidth=0.75)


def reset_coordinates(window, current_block):
    """
    If the unspanned block is connecting two things on different chromosomes,
    it will be necessary to set the start and ends to reflect this, as well as
    the orientation of the connection.
    """
    if current_block["region"]["start_chrom"] < window.chrom:
        current_block["region"]["start"] = 0
    if current_block["region"]["start_chrom"] > window.chrom:
        current_block["region"]["start"] = BIG_NUMBER

    if current_block["region"]["end_chrom"] > window.chrom:
        current_block["region"]["end"] = BIG_NUMBER
    if current_block["region"]["end_chrom"] < window.chrom:
        current_block["region"]["end"] = 0
    return current_block


def plot_unspanned_blocks(
    region_info_by_sample_order,
    max_sample_idx,
    block_count,
    spanned_block_ends,
    window,
    annotation_region_height,
    ax,
):
    """
    Go through all blocks in this regions and plot the ones that are not spanned.
    Use the stored spanned block ends to determine where the unspanned blocks start and end.
    """
    for sample_idx in range(max_sample_idx + 1):
        if sample_idx not in region_info_by_sample_order:
            continue
        matched_order_count = len(
            [
                r
                for r in region_info_by_sample_order[sample_idx]
                if len(r["coverages"]) == 0
            ]
        )
        if matched_order_count == 0:
            continue
        for i, breakend in enumerate(region_info_by_sample_order[sample_idx]):
            coverages = breakend["coverages"]
            if len(coverages) > 0:
                # spanned, already done
                continue
            sample_idx = breakend["sample_order_index"]

            plot_idx = (
                block_count + 1
            ) - sample_idx  # plot_idx is from top of image to bottom and shifted up to allow space for annotations
            if matched_order_count > 1:
                # subdivide plot index if multiple blocks in it
                plot_idx = (plot_idx - 1) + (i + 1) / matched_order_count

            plot_unspanned_block(
                breakend,
                plot_idx,
                window,
                spanned_block_ends,
                annotation_region_height,
                ax,
            )


def plot_spanned_blocks(
    region_info_by_sample_order,
    max_sample_idx,
    block_count,
    window,
    annotation_region_height,
    ax,
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
                block_count + 1
            ) - sample_idx  # plot_idx is from top of image to bottom and is shifted up by one to allow space for annotations
            if matched_order_count > 1:
                # subdivide plot index if multiple blocks in it
                plot_idx = (plot_idx - 1) + (i + 1) / matched_order_count

            # plots a block from the start of the alignment region to the end for a pair of
            # coordinates with spanning alignments, on the same level in the plot
            block_start_coords, block_end_coords = plot_spanned_block(
                coverages,
                plot_idx,
                window,
                orientation,
                matched_order_count,
                sample_idx,
                annotation_region_height,
                ax,
            )
            if len(block_start_coords) > 0:
                spanned_block_ends["starts"][tuple(block_start_coords)] = breakend
            if len(block_end_coords) > 0:
                spanned_block_ends["ends"][tuple(block_end_coords)] = breakend
    return spanned_block_ends


def plot_breaks(
    region_info_by_sample_order,
    max_sample_idx,
    window,
    block_count,
    window_count,
    annotation_region_height,
    ax,
):
    """
    Generates a static image representation of one complex SV.
    Starts with spanned blocks, then adds unspanned blocks.
    """
    spanned_block_ends = plot_spanned_blocks(
        region_info_by_sample_order,
        max_sample_idx,
        block_count,
        window,
        annotation_region_height,
        ax,
    )
    plot_unspanned_blocks(
        region_info_by_sample_order,
        max_sample_idx,
        block_count,
        spanned_block_ends,
        window,
        annotation_region_height,
        ax,
    )

    # configure ticks
    xtick_count = MAX_XTICK_COUNT
    for i in range(window_count):
        xtick_count -= 1
    padding = int(window.size / 10)
    ax.ticklabel_format(style="plain")
    ax.set_yticks([])
    rounding_factor = -3
    region_size = window.size + padding * 2
    if region_size < 1_000:
        rounding_factor = -2

    stepsize = max(1, int(region_size / xtick_count))
    xtick_locations = [
        round(x, rounding_factor)
        for x in range(window.start + padding, window.end - padding, stepsize)
    ]
    xtick_labels = [
        "{}".format(round_genomic_coords_to_str(0, int(l), max_at_kb=True))
        for l in xtick_locations
    ]
    ax.set_xticks(xtick_locations, xtick_labels, fontsize=FONTSIZE)

    window_length_annotation = round_genomic_coords_to_str(window.start, window.end)
    ax.set_title(
        "{}:{}-{} ({})".format(
            window.chrom,
            window.start,
            window.end,
            window_length_annotation,
        ),
        fontsize=FONTSIZE,
    )
    # add a little padding if there is none, to silence pyplot warnin
    if (window.start - padding) == (window.end + padding):
        padding += 1
    ax.set_xlim(window.start - padding, window.end + padding)
