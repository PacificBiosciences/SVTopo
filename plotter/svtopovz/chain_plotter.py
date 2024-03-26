#!/usr/bin/env python
from itertools import permutations

import sys
import matplotlib.pyplot as plt
import numpy as np
import logging

logger = logging.getLogger(__name__)

SMALL_BLOCK_SCALE = 30
ARROWHEAD_SCALE = 50
FONTSIZE = 7
MIN_PRECISION = 10


def plot_chain_representation(
    region_info_by_sample_order, coordinates, plot_height, axes
):
    """
    Plot the full structure of blocks in a chain format
    with ref versions and different sample versions for
    the possible paths through the region
    """
    region_size = coordinates["end"] - coordinates["start"]
    # MIN_PRECISION = region_size * MINIMUM_REGION_FRACT
    ref_blocks, ref_blocks_by_start, ref_blocks_by_end = generate_reference_blocks(
        region_info_by_sample_order
    )
    ref_blocks_by_idx = {i: b for b, i in ref_blocks.items()}
    skipped_ref_blocks = []
    for idx in ref_blocks_by_idx:
        ends = (ref_blocks_by_idx[idx][1], ref_blocks_by_idx[idx][3])
        if max(ends) - min(ends) < MIN_PRECISION:
            skipped_ref_blocks.append(idx)
    skipped_ref_blocks = sorted(skipped_ref_blocks)

    sample_paths = generate_sample_ordered_blocks(
        ref_blocks, ref_blocks_by_start, ref_blocks_by_end, region_info_by_sample_order
    )
    logger.debug("Sample structure paths: {}".format(sample_paths))
    if sample_paths is None:
        return False

    color_map = plt.get_cmap(name="tab20", lut=len(ref_blocks))
    colors = [color_map(i / len(ref_blocks)) for i in range(len(ref_blocks))]
    plot_reference_blocks(
        ref_blocks_by_idx,
        colors,
        plot_height,
        region_size,
        axes[0],
        skipped_ref_blocks,
    )
    return plot_sample_paths(
        sample_paths,
        ref_blocks_by_idx,
        colors,
        region_size,
        axes[1],
        skipped_ref_blocks,
    )


def generate_reference_blocks(region_info):
    """
    generate dict of blocks to letters with all regions represented in reference order
    """
    breaks = []
    for sample_order_idx in region_info:
        blocks = region_info[sample_order_idx]
        for block in blocks:
            breaks.append((block["region"]["start_chrom"], block["region"]["start"]))
            breaks.append((block["region"]["end_chrom"], block["region"]["end"]))
    breaks = sorted(list(set(breaks)), key=lambda b: (b[0], b[1]))

    expanded_blocks = {}
    expanded_blocks_by_start = {}
    expanded_blocks_by_end = {}
    ref_idx = 0
    for i in range(len(breaks) - 1):
        new_block = (breaks[i][0], breaks[i][1], breaks[i + 1][0], breaks[i + 1][1])
        expanded_blocks[new_block] = ref_idx
        ends = [(new_block[0], new_block[1]), (new_block[2], new_block[3])]
        ends = sorted(ends, key=lambda b: (b[0], b[1]))
        expanded_blocks_by_start[ends[0]] = new_block
        expanded_blocks_by_end[ends[1]] = new_block
        ref_idx += 1
    return expanded_blocks, expanded_blocks_by_start, expanded_blocks_by_end


def generate_sample_ordered_blocks(
    ref_blocks, ref_blocks_by_start, ref_blocks_by_end, region_info
):
    """
    Generate a list of lists, where each contains the reference numbers of the
    blocks in the sample and the orientations in the sample relative to the reference
    """
    max_tied_block_count = get_max_tied_block_count(region_info)
    if max_tied_block_count > 3:
        return
    sample_paths = [[] for _ in range(max_tied_block_count)]

    # previous connection is unspanned but unambiguous
    prev_info = {
        "is_direct_unspanned": True,
        "orientation": "",
        "start_index": None,
        "end_index": None,
    }
    for sample_order_idx in region_info:
        block_info = []
        for block in region_info[sample_order_idx]:
            block_region = block["region"]
            block_start = (block_region["start_chrom"], block_region["start"])
            block_end = (block_region["end_chrom"], block_region["end"])
            if block_start > block_end:
                tmp = block_end
                block_end = block_start
                block_start = tmp

            if block_start not in ref_blocks_by_start:
                continue
            block_start_idx = ref_blocks[ref_blocks_by_start[block_start]]

            if block_end not in ref_blocks_by_end:
                continue
            block_end_idx = ref_blocks[ref_blocks_by_end[block_end]]

            if len(block["coverages"]) > 0:
                new_spanned_blocks = get_spanned_blocks(
                    block, prev_info, block_start_idx, block_end_idx
                )
                block_info += new_spanned_blocks

            if len(block_info) > 0:
                prev_info["is_direct_unspanned"] = False
            else:
                prev_info["is_direct_unspanned"] = True

        prev_info["orientation"] = block["orientation"]
        prev_info["end_index"] = block_end_idx
        prev_info["start_index"] = block_start_idx

        if len(block_info) > 1:
            block_info = list(permutations(block_info))
        for i, sample_path in enumerate(sample_paths):
            if len(block_info) == 1:
                sample_path.append(block_info[0])
            elif len(block_info) == len(sample_paths):
                for block in block_info[i]:
                    sample_path.append(block)

    max_block_index = max(ref_blocks.values())
    extended_sample_paths = add_chain_ends(sample_paths, max_block_index)
    return extended_sample_paths


def get_max_tied_block_count(region_info):
    """
    Count the number of spanned blocks at each sample_order_index level and return the max
    """
    max_tied_block_count = 0
    for sample_order_idx in region_info:
        tied_block_count = 0
        for block in region_info[sample_order_idx]:
            if len(block["coverages"]) > 0:
                tied_block_count += 1
        if tied_block_count > max_tied_block_count:
            max_tied_block_count = tied_block_count
    return max_tied_block_count


def add_chain_ends(sample_paths: list, last_ref_block_index: int):
    """
    Add the ends of the chain so that it spans the entire region, not just the parts with
    rearrangements.
    """
    for sample_path in sample_paths:
        if len(sample_path[-1][0]) == 0:
            return

        # modify first block so that the chain starts at the beginning
        first_block_end = max(sample_path[0][0])
        new_first_block = (
            list(range(0, first_block_end + 1)),
            "+",
        )
        sample_path[0] = new_first_block

        # modify last block so that the chain goes all the way to the end
        last_block_start = min(sample_path[-1][0])
        new_last_block = (
            list(range(last_block_start, last_ref_block_index + 1)),
            "+",
        )
        sample_path[-1] = new_last_block
    return sample_paths


def get_spanned_blocks(
    block: dict, prev_info: dict, block_start_idx: int, block_end_idx: int
):
    """
    * `block` is the current block
    * `prev_info` is a dict with fields
        `is_direct_unspanned` (bool),
        `orientation` (string),
        `start_index` (int),
        `end_index` (int)
    * `block_start_idx` and `block_end_idx` are the indices of the start and end of this block within the chain
    """
    spanned_blocks = []
    if prev_info["is_direct_unspanned"]:
        if block["orientation"] == "":
            # if the current orientation is not known, it can be propagated from the previous block
            block["orientation"] = prev_info["orientation"]
        spanned_blocks.append(
            (
                list(range(block_start_idx, block_end_idx + 1)),
                block["orientation"],
            )
        )
    elif prev_info["end_index"] is not None:
        prev_leftmost = min(prev_info["start_index"], prev_info["end_index"])
        prev_rightmost = max(prev_info["start_index"], prev_info["end_index"])
        curr_leftmost = min(block_start_idx, block_end_idx)
        curr_rightmost = max(block_start_idx, block_end_idx)
        block_indices = []
        block_orientation = ""
        if prev_leftmost > curr_rightmost:
            # the previous block is to the right of the current one
            block_orientation = "-"
            block_indices = list(range(prev_leftmost - 1, curr_leftmost - 1, -1))

        else:
            # the previous block is to the left of the current one
            block_orientation = "+"
            block_indices = list(range(prev_rightmost, curr_rightmost + 1))

        for old_idx in (prev_info["start_index"], prev_info["end_index"]):
            if old_idx in block_indices:
                block_indices.remove(old_idx)
        spanned_blocks.append((block_indices, block_orientation))
    return spanned_blocks


def convert_to_letters(number, skipped_ref_blocks):
    """
    Converts the number indices to human-friendly letters, skipping blocks that are too small
    """
    if number < 0:
        sys.exit("Attempted to plot negative block")
    for skipped in skipped_ref_blocks:
        if skipped <= number:
            number -= 1
    result = ""
    while number >= 0:
        remainder = number % 26
        result = chr(ord("A") + remainder) + result
        number = number // 26 - 1
    return result


def plot_reference_blocks(
    ref_blocks,
    colors,
    plot_height,
    region_size,
    ax,
    skipped_ref_blocks,
):
    """
    plot genomic blocks in reference order
    """

    left_end = None
    for i in range(len(ref_blocks)):
        if i in skipped_ref_blocks:
            continue
        letter_code = convert_to_letters(i, skipped_ref_blocks)
        color = colors[i]
        patch_halfwidth = 0.25
        block_start = ref_blocks[i][1]
        if left_end is None or block_start < left_end:
            left_end = block_start
        block_end = ref_blocks[i][3]
        block_len = block_end - block_start
        if block_len < MIN_PRECISION:
            continue

        chain_middle = plot_height - patch_halfwidth

        if block_len < (region_size / SMALL_BLOCK_SCALE):
            # small block gets an arrowhead only
            block_shape = plt.Polygon(
                [
                    [block_start, chain_middle + patch_halfwidth],  # top left
                    [block_end, chain_middle],  # middle right
                    [block_start, chain_middle - patch_halfwidth],  # bottom left
                    [block_start, chain_middle + patch_halfwidth],  # top left
                ],
                fill=color,
                lw=0.5,
                color=color,
            )
            ax.add_patch(block_shape)
            ax.text(
                block_start + block_len / 2,
                chain_middle + patch_halfwidth * 2,
                letter_code,
                fontsize=FONTSIZE,
                horizontalalignment="right",
                verticalalignment="bottom",
            )
        else:
            arrowhead_len = region_size / ARROWHEAD_SCALE
            block_shape = plt.Polygon(
                [
                    [block_start, chain_middle + patch_halfwidth],  # top left
                    [
                        block_end - arrowhead_len,
                        chain_middle + patch_halfwidth,
                    ],  # arrowhead start
                    [block_end, chain_middle],  # arrowhead middle
                    [
                        block_end - arrowhead_len,
                        chain_middle - patch_halfwidth,
                    ],  # arrowhead end
                    [block_start, chain_middle - patch_halfwidth],  # bottom left
                    [block_start, chain_middle + patch_halfwidth],  # top left
                ],
                fill=color,
                lw=0.5,
                color=color,
            )
            ax.add_patch(block_shape)
            ax.text(
                block_start + ((block_len - arrowhead_len) / 2),
                chain_middle - (patch_halfwidth / 4),
                letter_code,
                fontsize=FONTSIZE,
                horizontalalignment="center",
                verticalalignment="center",
            )
    ax.text(
        left_end - (region_size / SMALL_BLOCK_SCALE) * 2,
        plot_height - 0.3,
        "Ref",
        fontsize=FONTSIZE,
    )


def plot_sample_paths(
    sample_paths,
    ref_blocks,
    colors,
    region_size,
    ax,
    skipped_ref_blocks,
):
    """
    Plots genomic blocks in each possible sample order
    """
    success = False
    max_option_plot_len = 0
    finished_paths = []
    # MIN_PRECISION = region_size * MINIMUM_REGION_FRACT
    for option_number, sample_path in enumerate(sample_paths):
        if sample_path in finished_paths:
            continue
        finished_paths.append(sample_path)
        option_plot_len = 0
        left_end = 0
        for block_numbers, orientation in sample_path:
            for block_number in block_numbers:
                letter_code = convert_to_letters(block_number, skipped_ref_blocks)
                color = colors[block_number]
                patch_halfwidth = 0.25
                block_start = min(
                    ref_blocks[block_number][1], ref_blocks[block_number][3]
                )
                block_end = max(
                    ref_blocks[block_number][1], ref_blocks[block_number][3]
                )
                block_len = block_end - block_start
                block_start = left_end
                block_end = block_start + block_len
                left_end = block_end
                option_plot_len += block_len
                chain_y_middle = len(sample_paths) - option_number + 0.5

                if orientation == "-":
                    polygon_coords = get_polygon(
                        block_len,
                        block_end,
                        block_start,
                        region_size,
                        chain_y_middle,
                        patch_halfwidth,
                        orientation,
                    )
                else:
                    polygon_coords = get_polygon(
                        block_len,
                        block_start,
                        block_end,
                        region_size,
                        chain_y_middle,
                        patch_halfwidth,
                        orientation,
                    )
                if block_len < MIN_PRECISION:
                    continue
                success = True
                if block_len < (region_size / SMALL_BLOCK_SCALE):
                    ax.text(
                        block_start + block_len / 2,
                        chain_y_middle + patch_halfwidth,
                        letter_code,
                        fontsize=FONTSIZE,
                        horizontalalignment="right",
                        verticalalignment="bottom",
                    )
                else:
                    arrowhead_len = region_size / ARROWHEAD_SCALE
                    plt.text(
                        block_start + ((block_len - arrowhead_len) / 2),
                        chain_y_middle - (patch_halfwidth / 4),
                        letter_code,
                        fontsize=FONTSIZE,
                        horizontalalignment="center",
                        verticalalignment="center",
                    )

                block_shape = plt.Polygon(
                    polygon_coords,
                    fill=color,
                    lw=0.5,
                    color=color,
                )
                plt.gca().add_patch(block_shape)
        if option_plot_len > max_option_plot_len:
            max_option_plot_len = option_plot_len
    if max_option_plot_len == 0:
        max_option_plot_len += 1

    if success:
        ax.set_ylim(0, 3)
        x_start = 0
        x_end = region_size
        if max_option_plot_len < region_size:
            difference = region_size - max_option_plot_len
            x_start = -difference
            x_end += region_size / difference
        elif max_option_plot_len > region_size:
            x_end = max_option_plot_len

        ax.set_xlim(x_start, x_end)
        if len(finished_paths) > 1:
            ax.set_title(
                "{} possible sample structures".format(len(finished_paths)),
                fontsize=FONTSIZE + 4,
            )
        else:
            ax.set_title("Sample structure", fontsize=FONTSIZE + 4)
    return success


def get_polygon(
    block_len,
    block_start,
    block_end,
    region_size,
    chain_middle,
    patch_halfwidth,
    orientation,
):
    """
    Get the coordinates to plot a polygon for a sample chain block
    """
    # default short block
    polygon_coords = [
        [block_start, chain_middle + patch_halfwidth],  # top left
        [block_end, chain_middle],  # middle right
        [
            block_start,
            chain_middle - patch_halfwidth,
        ],  # bottom left
        [block_start, chain_middle + patch_halfwidth],  # top left
    ]

    if block_len > (region_size / SMALL_BLOCK_SCALE):
        arrowhead_len = region_size / ARROWHEAD_SCALE
        if orientation == "-":
            arrowhead_len *= -1
        polygon_coords = [
            [block_start, chain_middle + patch_halfwidth],  # top left
            [
                block_end - arrowhead_len,
                chain_middle + patch_halfwidth,
            ],  # arrowhead start
            [block_end, chain_middle],  # arrowhead middle
            [
                block_end - arrowhead_len,
                chain_middle - patch_halfwidth,
            ],  # arrowhead end
            [
                block_start,
                chain_middle - patch_halfwidth,
            ],  # bottom left
            [block_start, chain_middle + patch_halfwidth],  # top left
        ]
    return polygon_coords
