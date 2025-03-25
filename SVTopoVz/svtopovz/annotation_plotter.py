import bisect
from collections import namedtuple
import matplotlib.pyplot as plt
import numpy as np
import logging

from svtopovz.utils import SMALL_FONTSIZE, BedRecord

logger = logging.getLogger(__name__)


def annotate(annotation_records, region, ax, level=0):
    """
    Add annotations to the plot for this region
    """
    assert level in [0, 1]
    region_length = region.end - region.start
    if region_length == 0:
        return

    if annotation_records is None:
        return
    region_bed_records = get_region_bed_records(annotation_records, region)
    if level == 0:
        annotation_top = 0.3
        annotation_bottom = 0.1
    elif level == 1:
        annotation_top = 1
        annotation_bottom = 0.8

    arrows_per_window = 100

    title_coords = []
    strand_coords = []
    for bed_record in region_bed_records:
        annotation_len = bed_record.end - bed_record.start
        annotation_fract = annotation_len / region_length

        add_box(bed_record, annotation_top, annotation_bottom, ax)
        add_strand(
            bed_record,
            annotation_len,
            annotation_fract,
            arrows_per_window,
            annotation_bottom,
            strand_coords,
            ax,
        )
        add_title(
            bed_record,
            annotation_len,
            annotation_top,
            annotation_bottom,
            title_coords,
            ax,
        )


def add_box(bed_record, annotation_top, annotation_bottom, ax):
    """
    Add the annotation of a bed_record rectangle
    """
    rectangle_corners = [
        [bed_record.start, annotation_top],
        [bed_record.start, annotation_bottom],
        [bed_record.end, annotation_bottom],
        [bed_record.end, annotation_top],
        [bed_record.start, annotation_top],
    ]
    rectangle = plt.Polygon(
        rectangle_corners,
        lw=0.25,
        color="lightblue",
    )
    ax.add_patch(rectangle)


def add_strand(
    bed_record,
    annotation_len,
    annotation_fract,
    arrows_per_window,
    annotation_bottom,
    strand_coords,
    ax,
):
    """
    Add the strand indicators to one bed_record rectangle
    """
    if annotation_fract == 0:
        return
    strand_indicator = {
        "": "",
        "+": ">",
        "-": "<",
    }
    arrow_count = int(annotation_fract * arrows_per_window)
    if arrow_count < 1:
        return
    for x_pos in np.arange(
        bed_record.start, bed_record.end, annotation_len / arrow_count
    ):
        overlaps_prev = False
        for strand_coord in strand_coords:
            if strand_coord[0] <= x_pos <= strand_coord[1]:
                overlaps_prev = True
        if not overlaps_prev:
            text = ax.text(
                x_pos,
                annotation_bottom,
                strand_indicator[bed_record.strand],
                fontsize=6,
            )
            plt.draw()

            # Get the bounding box of the text in display coordinates
            renderer = ax.figure.canvas.get_renderer()
            bbox = text.get_window_extent(renderer=renderer)

            # Convert the display coordinates back to data coordinates
            bbox_data = bbox.transformed(ax.transData.inverted())

            # Extract the x coordinates of the bounding box
            strand_coords.append((bbox_data.xmin, bbox_data.xmax))


def add_title(
    bed_record,
    annotation_len,
    annotation_top,
    annotation_bottom,
    title_coords,
    ax,
):
    """
    Add the title to one bed_record rectangle.
    Will truncate to not overgrow the size of the annotation rectangle.
    """
    annotation_x_middle = bed_record.start + (annotation_len / 2)
    annotation_title_y = annotation_bottom - (annotation_top - annotation_bottom)
    overlaps_prev = False
    for coord_pair in title_coords:
        if coord_pair[0] <= annotation_x_middle <= coord_pair[1]:
            overlaps_prev = True

    if not overlaps_prev:
        text = ax.text(
            annotation_x_middle,
            annotation_title_y,
            bed_record.title,
            fontsize=SMALL_FONTSIZE,
            horizontalalignment="center",
        )
        plt.draw()

        # Get the bounding box of the text in display coordinates
        renderer = ax.figure.canvas.get_renderer()
        bbox = text.get_window_extent(renderer=renderer)

        # Convert the display coordinates back to data coordinates
        bbox_data = bbox.transformed(ax.transData.inverted())

        # Extract the x coordinates of the bounding box
        title_coords.append((bbox_data.xmin, bbox_data.xmax))


def get_region_bed_records(bed_records, region):
    """
    If any of the bed records overlap the region, return them as list of
    truncated tuples (so the tuples lie only within the region).
    """
    region_records = []
    if bed_records is None:
        return
    chrom = (region.chrom).lower()
    if chrom not in bed_records:
        return region_records
    start = region.start
    end = region.end
    chrom_records = bed_records[chrom]
    region_start_idx = binary_get_start_index(chrom_records, start)

    for record in chrom_records[region_start_idx:]:
        if record.start > end:
            break
        new_record = BedRecord(
            max(start, record.start),
            min(end, record.end),
            record.strand,
            record.title,
        )
        region_records.append(new_record)
    return region_records


def binary_get_start_index(bed_records_list, start):
    """
    Binary search for the index of the start position or nearest lower value.
    Returns 0 if no lower value found
    """
    bisect_key_function = lambda x: x.end

    # find the index of the first item whose end is later than the start of the target region
    index = bisect.bisect_right(bed_records_list, start, key=bisect_key_function)

    if index is None:
        return 0
    return index
