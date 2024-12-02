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
import microhapdb
import pandas as pd
import sys
from util import load_markers, interval_distance


def main(markers, panel, ld_distance=10e6):
    final_panel_ids = extend_panel_iteratively(markers, panel.copy(), ld_distance=ld_distance)
    final_panel = markers[markers.Name.isin(final_panel_ids)]
    print_output(final_panel, panel)


def extend_panel_iteratively(markers, panel, ld_distance=10e6):
    extended_panel = panel.copy()
    while True:
        to_add = get_best_additions(markers, extended_panel, ld_distance=ld_distance)
        if len(to_add) == 0:
            break
        extended_panel |= to_add
    return extended_panel


def get_best_additions(markers, panel, ld_distance=10e6):
    to_add = set()
    candidates = get_candidates(markers, panel, ld_distance=ld_distance)
    ctable = markers[(markers.Name.isin(candidates)) & (markers.Ae >= 3.5)]
    for chrom, subtable in ctable.groupby("Chrom"):
        subtable = subtable.sort_values("Ae", ascending=False)
        to_add.add(subtable.Name.iloc[0])
    return to_add


def print_output(panel, original):
    print("Marker", "Chrom", "Location", "Extent", "Ae", "RSIDs", sep="\t")
    for mh in microhapdb.Marker.objectify(panel):
        print(mh.name, mh.chrom, mh.slug, len(mh), mh.data.Ae, ";".join(mh.varrefs), sep="\t")
    added = panel[~panel.Name.isin(original)][["Chrom", "Name", "Extent", "Ae"]]
    added.rename(columns={"Name": "Marker"})
    added.to_string(sys.stderr, index=False)
    print(f"\nAdded {len(added)} more microhaps for a total of {len(panel)} markers", file=sys.stderr)


def get_candidates(markers, panel, ld_distance=10e6):
    candidates = set()
    ingroup = markers[markers.Name.isin(panel)]
    outgroup = markers[~markers.Name.isin(panel)]
    for chrom, subingroup in ingroup.groupby("Chrom"):
        suboutgroup = outgroup[outgroup.Chrom == chrom]
        for i, outrow in suboutgroup.iterrows():
            for j, inrow in subingroup.iterrows():
                dist_ij = interval_distance((inrow.Start, inrow.End), (outrow.Start, outrow.End))
                if dist_ij < ld_distance:
                    break
            else:
                candidates.add(outrow.Name)
    return candidates


def get_parser():
    parser = ArgumentParser()
    parser.add_argument("markers", help="path to MicroHapDB marker definitions in CSV format")
    parser.add_argument("aes", help="path to MicroHapDB Ae table in CSV format")
    parser.add_argument("panel", help="path to a file with panel marker IDs")
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
    panel = pd.read_csv(args.panel, sep="\t")
    main(markers, set(panel.Marker), ld_distance=args.distance)
