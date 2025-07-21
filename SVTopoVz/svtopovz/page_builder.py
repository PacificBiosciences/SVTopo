#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create samplot vcf commands to execute and generate
companion HTML image browser.

Note: additional arguments are passed through to samplot plot
"""
from __future__ import print_function

import gzip
import logging
import os
import sys
from glob import glob
from .__init__ import __version__

from jinja2 import Environment, FileSystemLoader, select_autoescape
import shutil

try:
    from shlex import quote
except ImportError:
    from pipes import quote

logger = logging.getLogger(__name__)


def create_metadata(region_info):
    """
    creates a dict with the info about the SV
    that will be used in the website.
    """
    data_dict = {
        "chrom": region_info["chrom"],
        "chrom2": None,
        "start": region_info["pos"],
        "end": region_info["end"],
        "variant_ids": region_info["variant_ids"].sort(),
        "svlength": region_info["end"] - region_info["pos"],
        "samples": region_info["sample"],
        "nvariants": len(region_info["variant_ids"]),
        "image_name": region_info["image_name"],
    }
    return data_dict


def write_site(table_data, unique_table_data, out_dir, output_type):
    # grab the template
    template_dir = os.path.join(os.path.dirname(__file__), "templates")
    env = Environment(
        loader=FileSystemLoader(template_dir),
        autoescape=select_autoescape(["html"]),
    )

    # Add additional debug information if template cannot be found
    try:
        html_template = env.get_template("svtopo.html")
    except Exception as e:
        logger.error(f"Failed to load template: {e}")
        logger.error(f"Template directory: {template_dir}")
        logger.error(f"Directory exists: {os.path.exists(template_dir)}")
        if os.path.exists(template_dir):
            logger.error(f"Template directory contents: {os.listdir(template_dir)}")
        raise

    # Copy static artifacts to output directory
    for artifact_name in ["svtopo.png", "github-mark-white.png", "js/", "stylesheets/"]:
        src = os.path.join(template_dir, "static", artifact_name)
        dst = os.path.join(out_dir, artifact_name)
        if os.path.isdir(src):
            if os.path.isdir(dst):
                shutil.rmtree(dst)
            shutil.copytree(src, dst)
        else:
            shutil.copy2(src, dst)

    # write index.html
    with open("{out_dir}/index.html".format(out_dir=out_dir), "w") as fh:
        print(
            html_template.render(
                data=table_data,
                unique_data=unique_table_data,
                plot_type=output_type,
            ),
            file=fh,
        )


def readlines(file_path):
    """
    Generate lines from either a plain text file or gzipped
    """
    # Determine if the file is gzipped by its extension
    is_gzipped = file_path.endswith(".gz")

    # Open the file using gzip.open or built-in open based on the file type
    open_func = gzip.open if is_gzipped else open

    # Use context manager to handle file opening and closing
    with open_func(file_path, "rt") as f:
        for line in f:
            yield line.rstrip()  # Remove trailing newlines or spaces


def read_beds(svtopo_dir, image_type):
    # list of variant region info items
    image_info = []
    svtopo_beds = glob(os.path.join(svtopo_dir, "*bed"))
    svtopo_beds.extend(glob(os.path.join(svtopo_dir, "*bed.gz")))
    total_image_count = 0
    for bedfile in svtopo_beds:
        bed_prefix = ""
        finished_images = set()
        if bedfile.endswith(".gz"):
            bed_prefix = os.path.splitext(
                os.path.splitext(os.path.basename(bedfile))[0]
            )[0]
        else:
            bed_prefix = os.path.splitext(os.path.basename(bedfile))[0]
        is_svtopo_bed = False
        for line in readlines(bedfile):
            fields = line.strip().split()
            if line.startswith("#"):
                header_fields = line.strip().split()

                if len(fields) == 2 and header_fields[0] == "#SVTopo":
                    is_svtopo_bed = True
                    if header_fields[1] != "v{}".format(__version__):
                        logger.warning(
                            "SVTopo BED file {} version {} does not match current SVTopoVz version {}".format(
                                bedfile, header_fields[1], __version__
                            )
                        )
            elif is_svtopo_bed:
                if len(fields) < 4:
                    logger.error("SVTopo BED file {} is malformed".format(bedfile))
                    sys.exit()
                image_name = fields[3]
                if image_name in finished_images:
                    continue
                variant_ids = list(set(fields[4].split(","))) if len(fields) > 4 else []
                image_path = os.path.join(
                    svtopo_dir,
                    "images",
                    "{}_{}.{}".format(bed_prefix, image_name, image_type),
                )
                if not os.path.isfile(image_path):
                    continue
                for image_region in image_name.split("__"):
                    image_region_fields = image_region.split("-")
                    image_info_item = {
                        "sample": bed_prefix,
                        "chrom": image_region_fields[0],
                        "pos": int(image_region_fields[1]),
                        "end": int(image_region_fields[2]),
                        "image_name": "{}_{}".format(bed_prefix, image_name),
                        "variant_ids": variant_ids,
                    }
                    image_info.append(image_info_item)
                finished_images.add(image_name)
        logger.debug(
            "Generated {} images for prefix `{}`".format(
                len(finished_images), bed_prefix
            )
        )
        total_image_count += len(finished_images)
    logger.debug("Generated {} total images".format(total_image_count))
    return image_info


def generate_table(svtopo_dir, image_type):
    """
    Generate table data from bed files
    """
    table_data = []
    unique_table_data = []
    image_info = read_beds(svtopo_dir, image_type)
    included_images = set()
    for region in image_info:
        data_dict = create_metadata(region)
        table_data.append(data_dict)

        # store images by a single region
        image_name = region["image_name"]
        if image_name not in included_images:
            included_images.add(image_name)
            unique_table_data.append(data_dict)
    sort_key = lambda x: (
        x["chrom"] if x["chrom"] is not None else "",
        x["chrom2"] if x["chrom2"] is not None else "",
        x["start"] if x["start"] is not None else 0,
        x["end"] if x["end"] is not None else 0,
        x["samples"] if x["samples"] is not None else "",
        x["svlength"] if x["svlength"] is not None else 0,
        x["nvariants"] if x["nvariants"] is not None else 0,
        x["image_name"] if x["image_name"] is not None else "",
        str(x["variant_ids"]) if x["variant_ids"] is not None else "",
    )
    table_data.sort(key=sort_key)
    unique_table_data.sort(key=sort_key)

    return table_data, unique_table_data


def build_review_page(args):
    """
    Generate a review html page for the images created from one or more samples
    """
    table_data, unique_table_data = generate_table(args.svtopo_dir, args.image_type)
    write_site(table_data, unique_table_data, args.svtopo_dir, args.image_type)
    logger.info(
        "Open {}/index.html in your browser to see SVTopo Viewer".format(
            args.svtopo_dir.rstrip("/")
        )
    )
