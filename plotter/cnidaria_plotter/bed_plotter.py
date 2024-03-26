import bisect
from collections import namedtuple
import matplotlib.pyplot as plt
import numpy as np
import logging

logger = logging.getLogger(__name__)


def annotate(bed_records, region, ax):
    """
    Add bedfile annotations to the plot for this region
    """
    if bed_records is None:
        return
    region_bed_records = get_region_bed_records(bed_records, region)
    annotation_top = 0.4
    annotation_bottom = 0.2
    arrows_per_window = 50
    title_char_per_window = 175
    region_length = region["end"] - region["start"]
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
            ax,
        )
        add_title(
            bed_record,
            annotation_len,
            annotation_top,
            annotation_bottom,
            annotation_fract,
            title_char_per_window,
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
    ax,
):
    """
    Add the strand indicators to one bed_record rectangle
    """
    strand_indicator = {
        "": "",
        "+": ">",
        "-": "<",
    }

    arrow_count = int(annotation_fract * arrows_per_window)
    if arrow_count < 1:
        arrow_count = 1
    for x_pos in np.arange(
        bed_record.start, bed_record.end, annotation_len / arrow_count
    ):
        ax.text(
            x_pos,
            annotation_bottom,
            strand_indicator[bed_record.strand],
            fontsize=6,
        )


def add_title(
    bed_record,
    annotation_len,
    annotation_top,
    annotation_bottom,
    annotation_fract,
    chars_per_window,
    ax,
):
    """
    Add the title to one bed_record rectangle.
    Will truncate to not overgrow the size of the annotation rectangle.
    """
    annotation_x_middle = bed_record.start + (annotation_len / 2)
    annotation_title_y = annotation_bottom - (annotation_top - annotation_bottom)
    allowed_char_count = int(annotation_fract * chars_per_window)
    title = bed_record.title[:allowed_char_count]
    ax.text(
        annotation_x_middle,
        annotation_title_y,
        title,
        fontsize=4.5,
        horizontalalignment="center",
    )


def get_region_bed_records(bed_records, region):
    """
    If any of the bed records overlap the region, return them as list of
    truncated tuples (so the tuples lie only within the region).
    """
    region_records = []
    if bed_records is None:
        return
    chrom = (region["chrom"]).lower()
    if chrom not in bed_records:
        return region_records
    start = region["start"]
    end = region["end"]
    chrom_records = bed_records[chrom]
    region_start_idx = binary_get_start_index(chrom_records, start)

    bed_record = namedtuple("Record", ["start", "end", "strand", "title"])
    for record in chrom_records[region_start_idx:]:
        if record.start > end:
            break
        new_record = bed_record(
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
