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
import pandas as pd
import sys
from util import load_markers


def main(markers, panel, distance=10e6):
    panel_ids = panel.copy()
    while True:
        to_add = get_best_additions(markers, panel_ids, distance=distance)
        if len(to_add) == 0:
            break
        panel_ids |= to_add
    final_panel = markers[markers.Name.isin(panel_ids)][["Chrom", "Name", "Extent", "Ae"]]
    final_panel.to_csv(sys.stdout, sep=" ", index=False, header=False)
    added = final_panel[~final_panel.Name.isin(panel)]
    added.to_csv(sys.stderr, sep=" ", index=False, header=False)


def get_best_additions(markers, panel, distance=10e6):
    to_add = set()
    candidates = get_candidates(markers, panel, distance=distance)
    ctable = markers[(markers.Name.isin(candidates)) & (markers.Ae >= 3.5)]
    for chrom, subtable in ctable.groupby("Chrom"):
        if chrom == "chrX":
            continue
        subtable = subtable.sort_values("Ae", ascending=False)
        to_add.add(subtable.Name.iloc[0])
    return to_add


def get_candidates(markers, panel, distance=10e6):
    candidates = set()
    ingroup = markers[markers.Name.isin(panel)]
    outgroup = markers[~markers.Name.isin(panel)]
    for chrom, subingroup in ingroup.groupby("Chrom"):
        suboutgroup = outgroup[outgroup.Chrom == chrom]
        for i, outrow in suboutgroup.iterrows():
            smallest_distance = float("Inf")
            for j, inrow in subingroup.iterrows():
                x, y = sorted(((inrow.Start, inrow.End), (outrow.Start, outrow.End)))
                dist_xy = y[0] - x[1] if x[0] <= x[1] < y[0] else 0
                if dist_xy < smallest_distance:
                    smallest_distance = dist_xy
            if smallest_distance > distance:
                candidates.add(outrow.Name)
    return candidates


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
    main(markers, set(panel.Marker), distance=args.distance)
