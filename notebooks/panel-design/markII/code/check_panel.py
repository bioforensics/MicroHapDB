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
from itertools import combinations
import microhapdb
import pandas as pd
from util import load_markers, interval_distance


def main(markers, distance=10e6):
    for chrom, microhaps in markers.groupby("Chrom"):
        cmarkers = sorted(microhapdb.Marker.objectify(microhaps), key=lambda m: m.name)
        for mh1, mh2 in combinations(cmarkers, 2):
            mhdist = interval_distance((mh1.start, mh1.end), (mh2.start, mh2.end))
            if mhdist < distance:
                print(mh1.name, mh2.name, mhdist)


def get_parser():
    parser = ArgumentParser()
    parser.add_argument("markers", help="path to MicroHapDB marker definitions in CSV format")
    parser.add_argument("panel", help="path to a file with panel marker IDs")
    parser.add_argument("aes", help="path to MicroHapDB Ae table in CSV format")
    parser.add_argument(
        "--distance",
        type=float,
        default=10e6,
        metavar="D",
        help="minimum distance D between two markers to be considered independently inherited; by default D=10000000 (10 Mb)",
    )
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    markers = load_markers(args.markers, args.aes, objectify=False)
    panel = pd.read_csv(args.panel, sep=" ", names=["Chrom", "Marker", "Extent", "Ae"])
    panel_markers = markers[markers.Name.isin(panel.Marker)]
    main(panel_markers, distance=args.distance)
