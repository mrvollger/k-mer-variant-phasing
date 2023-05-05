#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Docs here
"""

__author__ = "Mitchell R. Vollger"
__credits__ = ["Mitchell R. Vollger"]
__maintainer__ = "Mitchell R. Vollger"
__email__ = "mrvollger@gmail.com"
__status__ = "Development"


import logging
import argparse
import pysam
import sys
import pandas as pd
from tqdm import tqdm


def parse():
    """Change all reads (and header) in a bam file to have one read group (RG)"""
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Input bam file",
        type=argparse.FileType("r"),
        default=sys.stdin,
    )
    parser.add_argument(
        "-r",
        "--read-list",
        help="List of read names and where to assign them",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--out",
        help="Output bam file.",
        default="-",
    )
    parser.add_argument(
        "-t", "--threads", help="Number of threads to use", type=int, default=8
    )
    parser.add_argument(
        "-v", "--verbose", help="increase logging verbosity", action="store_true"
    )
    args = parser.parse_args()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = logging.DEBUG if args.verbose else logging.WARNING
    logging.basicConfig(format=log_format, level=log_level)
    return args


def main():
    args = parse()
    bam = pysam.AlignmentFile(args.input, threads=args.threads, check_sq=False)
    out = sys.stdout if args.out == "-" else open(args.out, "w")
    out_bam = pysam.AlignmentFile(out, "wb", template=bam, threads=args.threads)
    reads = pd.read_csv(args.read_list, sep="\t")
    reads["HP"] = None
    reads["HP"]

    for rec in tqdm(bam.fetch(until_eof=True)):
        rec.set_tag("PS", 1)
        out_bam.write(rec)

    logging.info("Done")


if __name__ == "__main__":
    main()
