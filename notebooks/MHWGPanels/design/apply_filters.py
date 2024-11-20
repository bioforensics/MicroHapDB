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
from strategy import MaskingStrategy
import sys
from util import load_markers


def main(markers, masks, pass_stream=sys.stdout, fail_stream=None, vis_path=None):
    passing_markers = masks.apply_masks(markers).drop(columns=["Ae"])
    passing_markers.to_csv(pass_stream, index=False)
    if fail_stream:
        masks.masked_markers.to_csv(fail_stream, sep="\t", index=False)
    if vis_path:
        masks.visualize_masking(vis_path)


def get_parser():
    parser = ArgumentParser(description="Filter markers based on the provided mask(s)")
    parser.add_argument("markers", help="path to MicroHapDB marker definitions in CSV format")
    parser.add_argument("indel_mask", help="path to CSV file with indel mask")
    parser.add_argument("lowcomplex_mask", help="path to CSV file with low-complexity mask")
    parser.add_argument("repeat_mask", help="path to CSV file with repeat mask")
    parser.add_argument("str_mask", help="path to CSV file with forensic STR mask")
    parser.add_argument(
        "--aes",
        metavar="PATH",
        help="path to MicroHapDB Ae table in CSV format; by default, Ae values are reported as 0.0",
    )
    parser.add_argument(
        "--maxlen",
        type=int,
        metavar="L",
        default=260,
        help="filter out markers longer than L bp; by default L=260",
    )
    parser.add_argument(
        "--passed",
        metavar="FILE",
        default=sys.stdout,
        help="output file for markers passing filters; default is terminal (stdout)",
    )
    parser.add_argument("--failed", metavar="FILE", help="output file for markers failing filters")
    parser.add_argument(
        "--visualize",
        metavar="FILE",
        help="plot a graphic description of the masking results to FILE",
    )
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    markers = load_markers(args.markers, args.aes, objectify=False)
    masks = MaskingStrategy.from_tsv(
        args.indel_mask,
        args.lowcomplex_mask,
        args.repeat_mask,
        args.str_mask,
        max_length=args.maxlen,
    )
    main(markers, masks, pass_stream=args.passed, fail_stream=args.failed, vis_path=args.visualize)
