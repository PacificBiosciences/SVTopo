#!/usr/bin/env python
from itertools import permutations

import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from svtopovz.utils import *
import numpy as np
import logging

logger = logging.getLogger(__name__)


def plot_chain_representation(region_info_by_sample_order, windows, axes):
    """
    Plot the full structure of blocks in a chain format
    with ref versions and different sample versions for
    the possible paths through the region
    """
    ref_blocks, ref_blocks_by_start, ref_blocks_by_end = generate_reference_blocks(
        region_info_by_sample_order,
        windows,
    )
    total_region_size = 0
    for window in windows:
        total_region_size += window.size

    ref_blocks_by_idx = {i: b for b, i in ref_blocks.items()}
    skipped_ref_blocks = []
    for idx in ref_blocks_by_idx:
        ends = (ref_blocks_by_idx[idx][1], ref_blocks_by_idx[idx][3])
        if max(ends) - min(ends) < MIN_CLUSTER_PRECISION:
            skipped_ref_blocks.append(idx)
    skipped_ref_blocks = sorted(skipped_ref_blocks)

    sample_paths = generate_sample_ordered_blocks(
        ref_blocks, ref_blocks_by_start, ref_blocks_by_end, region_info_by_sample_order
    )
    logger.debug("Sample structure paths: {}".format(sample_paths))
    if sample_paths is None:
        logger.debug("No high-confidence paths found through sample structure")
        return False

    color_map = plt.get_cmap(name="tab20", lut=len(ref_blocks))
    colors = [color_map(i / len(ref_blocks)) for i in range(len(ref_blocks))]
    plot_reference_blocks_proportionally(
        ref_blocks_by_idx,
        colors,
        windows,
        axes[BREAKS],
        skipped_ref_blocks,
    )

    plot_reference_blocks_same_size(
        ref_blocks_by_idx,
        colors,
        total_region_size,
        axes[CHAINS],
        skipped_ref_blocks,
    )
    return plot_sample_paths(
        sample_paths,
        ref_blocks_by_idx,
        colors,
        total_region_size,
        axes[CHAINS],
        skipped_ref_blocks,
    )


def generate_reference_blocks(region_info, windows):
    """
    generate dict of blocks to letters with all regions represented in reference order
    """
    breaks = []
    for sample_order_idx in region_info:
        blocks = region_info[sample_order_idx]
        for block in blocks:
            breaks.append((block["region"]["start_chrom"], block["region"]["start"]))
            breaks.append((block["region"]["end_chrom"], block["region"]["end"]))
    window_breaks = []
    for window in windows:
        window_breaks.append((window.chrom, window.start))
        window_breaks.append((window.chrom, window.end))
    breaks += window_breaks[1:-1]
    breaks = sorted(list(set(breaks)), key=lambda b: (b[0], b[1]))

    expanded_blocks = {}
    expanded_blocks_by_start = {}
    expanded_blocks_by_end = {}
    ref_idx = 0
    for i in range(len(breaks) - 1):
        if breaks[i][0] != breaks[i + 1][0]:
            continue  # don't include the blocks that actually span between chromosomes
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
        "chrom": None,
    }
    for sample_order_idx in region_info:
        block_info = []
        block_found = False
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
            block_found = True

            if len(block["coverages"]) > 0:
                new_spanned_blocks = get_spanned_blocks(
                    block, prev_info, block_start_idx, block_end_idx
                )
                block_info += new_spanned_blocks

            if len(block_info) > 0:
                prev_info["is_direct_unspanned"] = False
            else:
                prev_info["is_direct_unspanned"] = True
        if block_found:
            prev_info["orientation"] = block["orientation"]
            prev_info["end_index"] = block_end_idx
            prev_info["start_index"] = block_start_idx
            prev_info["chrom"] = (
                block_region["start_chrom"]
                if block["orientation"] == "+"
                else block_region["end_chrom"]
            )

        if len(block_info) > 1:
            block_info = list(permutations(block_info))
        for i, sample_path in enumerate(sample_paths):
            if len(block_info) == 1:
                sample_path.append(block_info[0])
            elif len(block_info) == len(sample_paths):
                for block in block_info[i]:
                    sample_path.append(block)

    max_block_index = max(ref_blocks.values())
    extended_sample_paths = add_chain_ends(
        sample_paths, max_block_index, {v: k for k, v in ref_blocks.items()}
    )
    extended_sample_paths = flip_inverted_blocks(extended_sample_paths)
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


def add_chain_ends(
    sample_paths: list, last_ref_block_index: int, ref_blocks_by_idx: dict
):
    """
    Add the ends of the chain so that it spans the entire region, not just the parts with
    rearrangements.
    """
    for sample_path in sample_paths:
        if (
            len(sample_path) == 0
            or len(sample_path[0][0]) == 0
            or len(sample_path[-1][0]) == 0
        ):
            return

        # modify first block so that the chain starts at the beginning of the chromosome
        first_block_end_idx = max(sample_path[0][0])
        first_block_start_chrom = ref_blocks_by_idx[first_block_end_idx][0]
        new_first_block_idxs = []
        for i in range(0, first_block_end_idx + 1):
            if ref_blocks_by_idx[i][0] == first_block_start_chrom:
                new_first_block_idxs.append(i)

        new_first_block = (
            new_first_block_idxs,
            "+",
        )
        sample_path[0] = new_first_block

        # modify last block so that the chain goes all the way to the end of the chromosome
        last_block_start_idx = min(sample_path[-1][0])
        last_block_end_chrom = ref_blocks_by_idx[last_block_start_idx][0]
        new_last_block_idxs = []
        for i in range(last_block_start_idx, last_ref_block_index + 1):
            if ref_blocks_by_idx[i][0] == last_block_end_chrom:
                new_last_block_idxs.append(i)
        new_last_block = (
            new_last_block_idxs,
            "+",
        )
        sample_path[-1] = new_last_block
    return sample_paths

def flip_inverted_blocks(sample_paths):
    """
    Flips the ordering of sample indices for
    sample blocks that have a reverse orientation.
    """
    for i in range(len(sample_paths)):
        for j in range(len(sample_paths[i])):
            block_indices = sample_paths[i][j][0]
            block_orientation = sample_paths[i][j][1]
            if block_orientation == '-':
                block_indices = block_indices[::-1]
                block_info = [block_indices, block_orientation]
                sample_paths[i][j] = block_info
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
    if (
        prev_info["is_direct_unspanned"]
        or prev_info["chrom"] != block["region"]["start_chrom"]
    ):
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


def plot_reference_blocks_proportionally(
    ref_blocks,
    colors,
    windows,
    axes,
    skipped_ref_blocks,
):
    """
    Plot genomic blocks in reference order on the main axis,
    with sizing proportional to the actual block sizes
    """
    block_lens = []
    left_end = None
    for i, ax in enumerate(axes):
        window = windows[i]
        _, plot_height = ax.get_ylim()
        plot_height -= PATCH_HALFWIDTH * 1.5
        arrowhead_len = window.size / ARROWHEAD_SIZE
        if i == 0:
            ax.text(
                window.start - (arrowhead_len * 4),
                plot_height - 0.3,
                "Ref",
                fontsize=FONTSIZE,
            )
        for j, block_idx in enumerate(ref_blocks):
            if j in skipped_ref_blocks:
                continue
            letter_code = convert_to_letters(j, skipped_ref_blocks)
            color = colors[j]
            block = ref_blocks[block_idx]
            block_start = block[1]
            if left_end is None or block_start < left_end:
                left_end = block_start
            block_end = block[3]
            block_chrom = block[0]
            block_len = block_end - block_start
            if block_len < MIN_CLUSTER_PRECISION:
                continue
            starts_in_window = window.start <= block_start <= window.end
            ends_in_window = window.start <= block_end <= window.end
            block_lens.append((letter_code, block_len, color))
            if block_chrom != window.chrom or not (starts_in_window or ends_in_window):
                # if the block is outside this pane of the variant window, skip the lettering
                continue

            chain_middle = plot_height - PATCH_HALFWIDTH
            block_shape = plt.Polygon(
                [
                    [block_start, chain_middle + PATCH_HALFWIDTH],  # top left
                    [
                        block_end - arrowhead_len,
                        chain_middle + PATCH_HALFWIDTH,
                    ],  # arrowhead start
                    [block_end, chain_middle],  # arrowhead middle
                    [
                        block_end - arrowhead_len,
                        chain_middle - PATCH_HALFWIDTH,
                    ],  # arrowhead end
                    [block_start, chain_middle - PATCH_HALFWIDTH],  # bottom left
                    [block_start, chain_middle + PATCH_HALFWIDTH],  # top left
                ],
                fill=color,
                lw=0.5,
                color=color,
            )
            if block_len < (window.size / SMALL_BLOCK_SCALE):
                # small block gets an arrowhead only
                block_shape = plt.Polygon(
                    [
                        [block_start, chain_middle + PATCH_HALFWIDTH],  # top left
                        [block_end, chain_middle],  # middle right
                        [block_start, chain_middle - PATCH_HALFWIDTH],  # bottom left
                        [block_start, chain_middle + PATCH_HALFWIDTH],  # top left
                    ],
                    fill=color,
                    lw=0.5,
                    color=color,
                )
                ax.add_patch(block_shape)
                ax.text(
                    block_start + block_len / 2,
                    chain_middle + PATCH_HALFWIDTH,
                    letter_code,
                    fontsize=SMALL_FONTSIZE,
                    horizontalalignment="right",
                    verticalalignment="bottom",
                )
            else:
                ax.text(
                    block_start + ((block_len - arrowhead_len) / 2),
                    chain_middle - (PATCH_HALFWIDTH / 4),
                    letter_code,
                    fontsize=FONTSIZE,
                    horizontalalignment="center",
                    verticalalignment="center",
                )
            ax.add_patch(block_shape)

    plot_reference_block_sizes(block_lens)


def plot_reference_block_sizes(block_lens):
    """
    Given a list of tuples where each contains a block letter,length, and color,
    plot the sizes of the blocks on the right side of the main axis.
    """
    legend_elements = []
    finished = set()
    for block_letter, reference_block_len, color in block_lens:
        if block_letter in finished:
            continue
        finished.add(block_letter)
        block_len_str = round_genomic_coords_to_str(0, reference_block_len)
        label = "{}: {}".format(block_letter, block_len_str)
        legend_elements.append(mpatches.Patch(color=color, label=label))
    plt.legend(
        handles=legend_elements,
        loc="center left",
        bbox_to_anchor=(1, 3),
        fontsize=FONTSIZE,
    )


def plot_reference_blocks_same_size(
    ref_blocks,
    colors,
    region_size,
    ax,
    skipped_ref_blocks,
):
    """
    plot genomic blocks in reference order
    """
    block_len = min(region_size / len(ref_blocks), (region_size / (SMALL_BLOCK_SCALE)))

    left_end = 0
    for i in range(len(ref_blocks)):
        if i in skipped_ref_blocks:
            continue
        letter_code = convert_to_letters(i, skipped_ref_blocks)
        color = colors[i]
        chain_middle = 3.75  # top of the chain plot axis

        if i > 0 and ref_blocks[i - 1][0] != ref_blocks[i][0]:
            # add chromosomal transition
            ax.text(
                left_end,
                chain_middle - PATCH_HALFWIDTH,
                "//",
                fontsize=FONTSIZE + 4,
                # horizontalalignment="right",
                # verticalalignment="bottom",
            )
            left_end += block_len / 2.5

        if block_len < (region_size / SMALL_BLOCK_SCALE):
            # small block gets an arrowhead only
            block_shape = plt.Polygon(
                [
                    [left_end, chain_middle + PATCH_HALFWIDTH],  # top left
                    [left_end + block_len, chain_middle],  # middle right
                    [left_end, chain_middle - PATCH_HALFWIDTH],  # bottom left
                    [left_end, chain_middle + PATCH_HALFWIDTH],  # top left
                ],
                fill=color,
                lw=0.5,
                color=color,
            )
            ax.add_patch(block_shape)
            ax.text(
                left_end + block_len / 2,
                chain_middle + PATCH_HALFWIDTH * 2,
                letter_code,
                fontsize=FONTSIZE,
                horizontalalignment="right",
                verticalalignment="bottom",
            )
        else:
            arrowhead_len = region_size / ARROWHEAD_SIZE
            block_shape = plt.Polygon(
                [
                    [left_end, chain_middle + PATCH_HALFWIDTH],  # top left
                    [
                        (left_end + block_len) - arrowhead_len,
                        chain_middle + PATCH_HALFWIDTH,
                    ],  # arrowhead start
                    [(left_end + block_len), chain_middle],  # arrowhead middle
                    [
                        (left_end + block_len) - arrowhead_len,
                        chain_middle - PATCH_HALFWIDTH,
                    ],  # arrowhead end
                    [left_end, chain_middle - PATCH_HALFWIDTH],  # bottom left
                    [left_end, chain_middle + PATCH_HALFWIDTH],  # top left
                ],
                fill=color,
                lw=0.5,
                color=color,
            )
            ax.add_patch(block_shape)
            ax.text(
                left_end + ((block_len - arrowhead_len) / 2),
                chain_middle - (PATCH_HALFWIDTH / 4),
                letter_code,
                fontsize=FONTSIZE,
                horizontalalignment="center",
                verticalalignment="center",
            )
        left_end += block_len
    ax.text(
        0,
        chain_middle + (PATCH_HALFWIDTH * 2),
        "Reference path:",
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
    Plots genomic blocks in each possible sample order up to a maximum of three options
    """
    successes = set()
    max_option_plot_len = 0
    finished_paths = []
    for path_number, sample_path in enumerate(sample_paths):
        if sample_path in finished_paths:
            continue
        if len(sample_paths) > 1:
            chain_y_middle = len(sample_paths) - (path_number * 1.25)
        else:
            chain_y_middle = 1.5

        finished_paths.append(sample_path)
        option_plot_len = 0
        left_end = 0
        for block_numbers, orientation in sample_path:
            for i, block_number in enumerate(block_numbers):
                if block_number in skipped_ref_blocks:
                    continue
                option_plot_len, success, left_end = plot_block(
                    ax,
                    block_numbers,
                    i,
                    ref_blocks,
                    skipped_ref_blocks,
                    colors,
                    chain_y_middle,
                    region_size,
                    orientation,
                    option_plot_len,
                    left_end,
                )
                successes.add(success)

        label = "Sample path {}:"
        ax.text(
            0,
            chain_y_middle + (PATCH_HALFWIDTH * 2),
            label.format(path_number),
            fontsize=FONTSIZE,
        )
        if option_plot_len > max_option_plot_len:
            max_option_plot_len = option_plot_len
    if max_option_plot_len == 0:
        max_option_plot_len += 1

    if True in successes:
        ax.set_ylim(0, 4)
        x_end = max(region_size, max_option_plot_len)
        ax.set_xlim(0, x_end)
    return success


def plot_block(
    ax,
    block_numbers,
    path_index,
    ref_blocks,
    skipped_ref_blocks,
    colors,
    chain_y_middle,
    region_size,
    orientation,
    option_plot_len,
    left_end,
):
    """
    Plot a block in a sample path through the variant graph
    """
    block_len = min(region_size / len(ref_blocks), (region_size / (SMALL_BLOCK_SCALE)))
    block_number = block_numbers[path_index]

    if (
        path_index > 0
        and block_numbers[path_index - 1] == (block_numbers[path_index] - 1)
        and ref_blocks[block_numbers[path_index - 1]][0]
        != ref_blocks[block_numbers[path_index]][0]
    ):
        # add chromosomal transition
        ax.text(
            left_end,
            chain_y_middle - PATCH_HALFWIDTH,
            "//",
            fontsize=FONTSIZE + 4,
        )
        left_end += block_len / 2.5

    letter_code = convert_to_letters(block_number, skipped_ref_blocks)
    color = colors[block_number]
    block_start = min(ref_blocks[block_number][1], ref_blocks[block_number][3])
    block_end = max(ref_blocks[block_number][1], ref_blocks[block_number][3])
    block_start = left_end
    block_end = block_start + block_len

    left_end = block_end
    option_plot_len += block_len
    arrowhead_len = region_size / ARROWHEAD_SIZE

    if orientation == "-":
        arrowhead_len *= -1
        polygon_coords = get_polygon(
            block_len,
            block_end,
            block_start,
            region_size,
            chain_y_middle,
            orientation,
        )
    else:
        polygon_coords = get_polygon(
            block_len,
            block_start,
            block_end,
            region_size,
            chain_y_middle,
            orientation,
        )
    if block_len < MIN_CLUSTER_PRECISION:
        return option_plot_len, False, left_end
    if block_len >= (region_size / SMALL_BLOCK_SCALE):
        ax.text(
            block_start + ((block_len - arrowhead_len) / 2),
            chain_y_middle - (PATCH_HALFWIDTH / 4),
            letter_code,
            fontsize=FONTSIZE,
            horizontalalignment="center",
            verticalalignment="center",
        )
    else:
        ax.text(
            block_start + block_len / 2,
            chain_y_middle + PATCH_HALFWIDTH,
            letter_code,
            fontsize=FONTSIZE,
            horizontalalignment="right",
            verticalalignment="bottom",
        )

    block_shape = plt.Polygon(
        polygon_coords,
        fill=color,
        lw=0.5,
        color=color,
    )
    ax.add_patch(block_shape)
    return option_plot_len, True, left_end


def get_polygon(
    block_len,
    block_start,
    block_end,
    region_size,
    chain_middle,
    orientation,
):
    """
    Get the coordinates to plot a polygon for a sample chain block
    """
    # default short block
    polygon_coords = [
        [block_start, chain_middle + PATCH_HALFWIDTH],  # top left
        [block_end, chain_middle],  # middle right
        [
            block_start,
            chain_middle - PATCH_HALFWIDTH,
        ],  # bottom left
        [block_start, chain_middle + PATCH_HALFWIDTH],  # top left
    ]

    if block_len >= (region_size / SMALL_BLOCK_SCALE):
        arrowhead_len = region_size / ARROWHEAD_SIZE
        if orientation == "-":
            arrowhead_len *= -1
        polygon_coords = [
            [block_start, chain_middle + PATCH_HALFWIDTH],  # top left
            [
                block_end - arrowhead_len,
                chain_middle + PATCH_HALFWIDTH,
            ],  # arrowhead start
            [block_end, chain_middle],  # arrowhead middle
            [
                block_end - arrowhead_len,
                chain_middle - PATCH_HALFWIDTH,
            ],  # arrowhead end
            [
                block_start,
                chain_middle - PATCH_HALFWIDTH,
            ],  # bottom left
            [block_start, chain_middle + PATCH_HALFWIDTH],  # top left
        ]
    return polygon_coords
