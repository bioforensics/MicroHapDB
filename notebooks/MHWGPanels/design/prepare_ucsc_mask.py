#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd
import sys


def main(instream, outstream, sine=929, line=411, ltr=909):
    repeats = parse_ucsc_rmsk_track(instream)
    repeats = repeats[
        ((repeats.repClass == "SINE") & (repeats.swScore > sine))
        | ((repeats.repClass == "LINE") & (repeats.swScore > line))
        | ((repeats.repClass == "LTR") & (repeats.swScore > ltr))
    ]
    repeats = repeats[["genoName", "genoStart", "genoEnd"]]
    repeats.to_csv(outstream, sep="\t", index=False, header=False)


def parse_ucsc_rmsk_track(path):
    header = [
        "bin",
        "swScore",
        "milliDiv",
        "milliDel",
        "milliIns",
        "genoName",
        "genoStart",
        "genoEnd",
        "genoLeft",
        "strand",
        "repName",
        "repClass",
        "repFamily",
        "repStart",
        "repEnd",
        "repLeft",
        "id",
    ]
    return pd.read_csv(path, sep="\t", names=header)


def get_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "-o",
        "--out",
        metavar="PATH",
        help="write output to PATH; by default, output is written to terminal (stdout)",
    )
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
    parser.add_argument("rmsk", help="path to RepeatMasker annotation file")
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    if args.out is None:
        args.out = sys.stdout
    main(args.rmsk, args.out, sine=args.sine, line=args.line, ltr=args.ltr)
