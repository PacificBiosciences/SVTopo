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
    if not path.exists(filepath) and path.isfile(filepath):
        return False
    return filepath


def valid_prefix(prefix):
    if not path.exists(path.dirname(prefix)) and path.isfile(path.dirname(prefix)):
        return False
    return prefix


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
        "-j",
        "--json",
        help="path to json-formatted complex SV data, extracted from BAM file using SVTopo. GZIP allowed.",
        required=True,
        type=valid_file,
    )
    parser.add_argument(
        "-p",
        "--out-prefix",
        help="prefix for output files. May include a directory path.",
        required=valid_prefix,
    )
    parser.add_argument(
        "-b",
        "--bed",
        help="paths to genome annotations in BED file format",
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
        help="maximum gap size to show in one panel, in megabases. Default is entire chromosome",
        default=0,
        type=float,
    )
    parser.add_argument(
        "--log-level",
        help="set log level",
        choices=["debug", "info"],
        default="info",
    )
    parser.add_argument(
        "--ignore-simple-dels",
        help="identifies simple deletions (defined as having exactly two forward spanned blocks, one forward unspanned) and skips them",
        action="store_true",
    )
    parser.add_argument(
        "--ignore-simple-dups",
        help="identifies simple duplications (defined as having exactly two forward spanned blocks, one reverse unspanned) and skips them",
        action="store_true",
    )
    return parser


def main():
    print("\nSVTopoVz v{}".format(__version__), file=sys.stderr)
    parser = setup_args()
    args = parser.parse_args()

    if args.log_level == "info":
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.DEBUG)

    prefix_path = args.out_prefix
    prefix_dir = path.dirname(prefix_path)
    if len(prefix_dir) > 0 and not path.isdir(prefix_dir):
        # there is a directory included in the prefix path but it doesn't exist
        logger.error(
            " prefix includes directory `{}`, which does not exist".format(prefix_dir)
        )
        sys.exit()

    svtopovz(args)


if __name__ == "__main__":
    print(
        "You are running this module directly, which should only be done for debugging",
        file=sys.stderr,
    )
    sys.exit(main() or 0)
