#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd
import re
import sys


def main(instream, outstream, distance=10e6):
    strs = parse_fssg_info(instream)
    strs["Start"] = strs.Start.apply(lambda x: max(0, x - int(distance)))
    strs["End"] = strs.End + int(distance)
    strs = strs[["Chrom", "Start", "End"]]
    strs.to_csv(outstream, sep="\t", index=False, header=False)


def parse_fssg_info(path):
    table = pd.read_csv(path)
    table["Chrom"] = table.FullRangeGRCh38.apply(lambda x: parse_range(x)[0])
    table["Start"] = table.FullRangeGRCh38.apply(lambda x: parse_range(x)[1])
    table["End"] = table.FullRangeGRCh38.apply(lambda x: parse_range(x)[2])
    return table


def parse_range(chrom_range):
    match = re.match(r"(chr\d+):(\d+)-(\d+)", chrom_range)
    chrom, start, end = match.groups()
    return chrom, int(start), int(end)


def get_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "-o",
        "--out",
        metavar="PATH",
        help="write output to PATH; by default, output is written to terminal (stdout)",
    )
    parser.add_argument(
        "--distance",
        type=float,
        default=10e6,
        metavar="D",
        help="extend beyond the end of each locus D bp; by default D=10000000",
    )
    parser.add_argument("fssg", help="path to FSSG info table")
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    if args.out is None:
        args.out = sys.stdout
    main(args.fssg, args.out, distance=args.distance)
