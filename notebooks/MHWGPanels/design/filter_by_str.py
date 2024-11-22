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
import re
from tqdm import tqdm
from util import load_markers


def main(markers, strs):
    print("Marker", "Extent", "Ae", "Locus", sep="\t")
    for mh in (pbar := tqdm(markers)):
        pbar.set_description(f"{mh.name:<20}")
        overlap = strs[(strs.Chrom == mh.chrom) & (strs.End > mh.start) & (strs.Start < mh.end)]
        if len(overlap) > 0:
            loci = ";".join(sorted(overlap.Locus))
            print(mh.name, len(mh), mh.data.Ae, loci, sep="\t")


def parse_fssg_info(path, distance=10e6):
    table = pd.read_csv(path)
    table["Chrom"] = table.FullRangeGRCh38.apply(lambda x: parse_range(x)[0])
    table["Start"] = table.FullRangeGRCh38.apply(lambda x: parse_range(x)[1])
    table["End"] = table.FullRangeGRCh38.apply(lambda x: parse_range(x)[2])
    table["Start"] = table.Start.apply(lambda x: max(0, x - int(distance)))
    table["End"] = table.End + int(distance)
    return table


def parse_range(chrom_range):
    match = re.match(r"(chr\d+):(\d+)-(\d+)", chrom_range)
    chrom, start, end = match.groups()
    return chrom, int(start), int(end)


def get_parser():
    parser = ArgumentParser()
    parser.add_argument("markers", help="path to MicroHapDB marker definitions in CSV format")
    parser.add_argument("fssg", help="path to FSSG info table")
    parser.add_argument(
        "--distance",
        type=float,
        default=10e6,
        metavar="D",
        help="any microhap within D bp of an FSSG STR marker will be considered in linkage disequilibrium with that STR and filtered out; by default D=10000000 (10 Mb)",
    )
    parser.add_argument(
        "--aes",
        metavar="PATH",
        help="path to MicroHapDB Ae table in CSV format; by default, Ae values are reported as 0.0",
    )
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    markers = load_markers(args.markers, args.aes)
    strs = parse_fssg_info(args.fssg, distance=args.distance)
    main(markers, strs)
