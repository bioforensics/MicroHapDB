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
import polars as pl
from tqdm import tqdm
from util import load_markers, parse_ucsc_rmsk_track


def main(markers, repeats, distance=4):
    print("Marker", "Extent", "Ae", "Offset", "Repeats", sep="\t")
    for mh in tqdm(markers):
        offending_offsets = set()
        offending_repeats = set()
        for offset in mh.offsets:
            start = offset - distance
            end = offset + distance
            overlapping = repeats.filter(
                (pl.col("genoName") == mh.chrom)
                & (pl.col("genoStart") < end)
                & (start < pl.col("genoEnd"))
            )
            if len(overlapping) > 0:
                offending_offsets.add(offset)
                offending_repeats.update(overlapping.get_column("repName"))
        if len(offending_offsets) > 0:
            offset_data = ";".join(sorted(map(str, offending_offsets)))
            repeat_data = ";".join(sorted(offending_repeats))
            print(mh.name, len(mh), mh.data.Ae, offset_data, repeat_data, sep="\t")


def get_parser():
    parser = ArgumentParser()
    parser.add_argument("markers", help="path to MicroHapDB marker definitions in CSV format")
    parser.add_argument("aes", help="path to MicroHapDB Ae table in CSV format")
    parser.add_argument("rmsk", help="path to UCSC RepeatMasker track data")
    parser.add_argument(
        "--distance",
        type=int,
        default=4,
        metavar="D",
        help="flag markers with ADSs within D bp of a dbSNP indel; by default D=4",
    )
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    markers = load_markers(args.markers, args.aes)
    repeats = parse_ucsc_rmsk_track(args.rmsk)
    repeats = repeats.filter(pl.col("repClass") == "Low_complexity")
    main(markers, repeats, distance=args.distance)
