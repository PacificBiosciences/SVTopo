#!/usr/bin/env python
from collections import namedtuple

import gzip
import json
import os
import sys
from cnidaria_plotter import bnd_plotter
import logging

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
            return json.load(json_fh)
    else:
        with open(json_filename, "r") as json_fh:
            return json.load(json_fh)


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


def cnidaria_plotter(args):
    sv_info = unpack_json(args.json)
    bed_records = unpack_bed_records(args.bed)
    for count, event_info in enumerate(sv_info):
        bnd_plotter.plot_region(
            event_info, args.out_prefix, bed_records, args.image_type
        )
        logger.debug("Finished event #{}".format(count))
