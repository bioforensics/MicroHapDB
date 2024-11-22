#!/usr/bin/env python

# -------------------------------------------------------------------------------------------------
# Copyright (c) 2024, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from argparse import ArgumentParser
from linkage_graph import LinkageGraph, logger
import sys
from util import load_markers


def main(markers, max_per_chrom_long=8, max_per_chrom_short=8, batches=8, distance=10e6):
    logger.info(f"Populating linkage graph long={max_per_chrom_long} short={max_per_chrom_short}")
    graph = LinkageGraph.populate(
        markers,
        max_per_chrom_long=max_per_chrom_long,
        max_per_chrom_short=max_per_chrom_short,
        batches=batches,
        distance=distance,
    )
    panel = graph.panel
    for mh in panel:
        print(mh.chrom, mh.name, len(mh), mh.data.Ae)
    print(len(panel), sep="\n", file=sys.stderr)


def get_parser():
    parser = ArgumentParser()
    parser.add_argument("markers", help="path to MicroHapDB marker definitions in CSV format")
    parser.add_argument("aes", help="path to MicroHapDB Ae table in CSV format")
    parser.add_argument(
        "--max-per-chrom",
        type=int,
        nargs=2,
        default=(8, 8),
        metavar="M",
        help="select the M highest ranked short and long markers (respectively) per chromosome by Ae; by default M=(8, 8)",
    )
    parser.add_argument(
        "--batches",
        type=int,
        default=8,
        metavar="B",
        help="split the linkage graph into B batches; by default B=8",
    )
    parser.add_argument(
        "--distance",
        type=float,
        default=10e6,
        metavar="D",
        help="two markers must be separated by more than D bp to be considered independently inherited (as a heuristic); by default D=10000000",
    )
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    markers = load_markers(args.markers, args.aes)
    main(
        markers,
        max_per_chrom_short=args.max_per_chrom[0],
        max_per_chrom_long=args.max_per_chrom[1],
        batches=args.batches,
        distance=args.distance,
    )
