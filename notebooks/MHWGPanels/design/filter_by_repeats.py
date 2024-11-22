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
import sys
from tqdm import tqdm
from util import load_markers, parse_ucsc_rmsk_track


def main(markers, repeats, whitelist):
    print("Marker", "Extent", "Ae", "Repeats", sep="\t")
    for mh in (pbar := tqdm(markers)):
        pbar.set_description(f"{mh.name:<20}")
        if whitelist and mh.locus in whitelist:
            continue
        overlap = repeats.filter(
            (pl.col("genoName") == mh.chrom)
            & (pl.col("genoEnd") > mh.start)
            & (pl.col("genoStart") < mh.end)
        )
        if len(overlap) > 0:
            print(mh.name, len(mh), mh.data.Ae, format_repeat_results(overlap), sep="\t")


def format_repeat_results(results):
    slugs = list()
    for repeat in results.iter_rows(named=True):
        slug = f"{repeat['repClass']}>{repeat['repFamily']}>{repeat['repName']}[{repeat['genoName']}:{repeat['genoStart']}-{repeat['genoEnd']}]"
        slugs.append(slug)
    return ";".join(slugs)


def load_and_filter_repeats(rmsk_path, sine=929, line=411, ltr=909):
    repeats = parse_ucsc_rmsk_track(rmsk_path)
    print(f"Loaded RepeatMasker data; {len(repeats)} annotated repeats", file=sys.stderr)
    repeats = repeats.filter(
        ((pl.col("repClass") == "SINE") & (pl.col("swScore") > sine))
        | ((pl.col("repClass") == "LINE") & (pl.col("swScore") > line))
        | ((pl.col("repClass") == "LTR") & (pl.col("swScore") > ltr))
    )
    print(f"Filtered repeats; {len(repeats)} repeats remain", file=sys.stderr)
    return repeats


def get_parser():
    parser = ArgumentParser()
    parser.add_argument("markers", help="path to MicroHapDB marker definitions in CSV format")
    parser.add_argument("rmsk", help="path to UCSC RepeatMasker track data")
    parser.add_argument(
        "--sine",
        type=float,
        default=929,
        metavar="X",
        help="discard SINE elements with a RepeatMasker score of X or lower; by default X=929",
    )
    parser.add_argument(
        "--line",
        type=float,
        default=411,
        metavar="Y",
        help="discard LINE elements with a RepeatMasker score of Y or lower; by default Y=411",
    )
    parser.add_argument(
        "--ltr",
        type=float,
        default=909,
        metavar="Z",
        help="discard LTR elements with a RepeatMasker score of Z or lower; by default Z=909",
    )
    parser.add_argument(
        "--aes",
        metavar="FILE",
        help="path to MicroHapDB Ae table in CSV format; by default, Ae values are reported as 0.0",
    )
    parser.add_argument(
        "--whitelist",
        metavar="FILE",
        help="path to a text file with locus IDs indicating loci for which primers have been successfully designed previously",
    )
    return parser


def load_whitelist(path):
    with open(path, "r") as fh:
        return set(fh.read().split())


if __name__ == "__main__":
    args = get_parser().parse_args()
    markers = load_markers(args.markers, args.aes)
    repeats = load_and_filter_repeats(args.rmsk, sine=args.sine, line=args.line, ltr=args.ltr)
    whitelist = None if args.whitelist is None else load_whitelist(args.whitelist)
    main(markers, repeats, whitelist)
