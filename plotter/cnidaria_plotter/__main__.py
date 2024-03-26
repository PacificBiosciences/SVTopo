#!/usr/bin/env python
from __future__ import print_function

import argparse
import sys
from os import path
import logging

logger = logging.getLogger(__name__)

from .__init__ import __version__
from .cnidaria_plotter import cnidaria_plotter


def setup_args():
    parser = argparse.ArgumentParser(
        prog="cnidaria_plotter", formatter_class=argparse.ArgumentDefaultsHelpFormatter
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
        help="path to json-formatted complex SV data, extracted from BAM file using cnidaria. GZIP allowed.",
        required=True,
    )
    parser.add_argument(
        "-p",
        "--out-prefix",
        help="prefix for output files. May include a directory path.",
        required=True,
    )
    parser.add_argument(
        "-b",
        "--bed",
        help="paths to genome annotations in BED file format",
        required=False,
    )
    parser.add_argument(
        "--image-type",
        help="type of image to generate",
        default="png",
        choices=["png", "jpg", "jpeg", "svg", "pdf"],
    )
    parser.add_argument(
        "--log-level",
        help="set log level",
        choices=["debug", "info"],
        default="info",
    )
    return parser


def main():
    print("\nCNIDARIA_PLOTTER v{}".format(__version__), file=sys.stderr)
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

    cnidaria_plotter(args)


if __name__ == "__main__":
    print(
        "You are running this module directly, which should only be done for debugging",
        file=sys.stderr,
    )
    sys.exit(main() or 0)
