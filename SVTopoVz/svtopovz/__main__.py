#!/usr/bin/env python
from __future__ import print_function

import argparse
import sys
from os import path
import logging

logger = logging.getLogger(__name__)

from .__init__ import __version__
from .svtopovz import svtopovz


def valid_file(filepath):
    if not path.isfile(filepath):
        raise argparse.ArgumentTypeError(f"{filepath} is not a valid file path.")
    return filepath


def valid_dir(dirpath):
    if not path.isdir(dirpath):
        raise argparse.ArgumentTypeError(f"{dirpath} is not a valid directory path.")
    return dirpath


def positive_int(value):
    """
    Defines an integer as positive
    """
    value = int(value)
    if value <= 0:
        raise argparse.ArgumentTypeError(f"{value} is not a valid positive integer.")
    else:
        return value


def setup_args():
    parser = argparse.ArgumentParser(
        prog="svtopovz", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-v",
        "--version",
        help="Installed version ({})".format(__version__),
        action="version",
        version="%(prog)s " + str(__version__),
    )
    parser.add_argument(
        "--svtopo-dir",
        help="path to directory containing one or more svtopo output file pairs (in json and bed format). GZIP allowed. ",
        required=True,
        type=valid_dir,
    )
    parser.add_argument(
        "--annotation-bed",
        help="space delimited list of one or more paths to genome annotations in BED file format - optionally allows annotation title in column 4 and strand (+/-) in column 5",
        required=False,
        type=valid_file,
        nargs="*",
    )
    parser.add_argument(
        "--genes",
        help="single path to gene annotations in GFF3 or GTF format (based on GENCODE v45 annotations)",
        required=False,
        type=valid_file,
    )
    parser.add_argument(
        "--image-type",
        help="type of image to generate",
        default="png",
        choices=["png", "jpg", "jpeg", "svg", "pdf"],
    )
    parser.add_argument(
        "--max-gap-size-mb",
        help="maximum gap size to show in one panel, in megabases. Default is 0.5mB",
        default=0.5,
        type=float,
    )
    parser.add_argument(
        "--verbose",
        help="print verbose output for debugging purposes",
        action="store_true",
    )
    parser.add_argument(
        "--include-simple-breakpoints",
        help="does not skip simple single-breakpoint events, such as deletions, duplications, and nonreciprocal translocations",
        action="store_true",
    )
    return parser


def main():
    print("\nSVTopoVz v{}".format(__version__), file=sys.stderr)
    parser = setup_args()
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if not path.isdir(args.svtopo_dir):
        logger.error("SVTopo directory {} does not exist".format(args.svtopo_dir))
        sys.exit()

    svtopovz(args)


if __name__ == "__main__":
    print(
        "You are running this module directly, which should only be done for debugging",
        file=sys.stderr,
    )
    sys.exit(main() or 0)
