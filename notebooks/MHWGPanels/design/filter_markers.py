#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd


def get_parser():
    parser = ArgumentParser(description="Discard markers based on the provided mask(s)")
    parser.add_argument("markers", help="MicroHapDB markers in CSV format")
    parser.add_argument("-o", "--out", metavar="FILE", help="output CSV file containing markers passing filters; by default output is written to the terminal (stdout)")
    parser.add_argument("-m", "--marker-mask", metavar="BED", nargs="+", help="path(s) to one or more BED files; markers overlapping with the specified intervals will be discarded")
    parser.add_argument("-s", "--snp-mask", metavar="BED", nargs="+", help="path(s) to one or more BED files; markers with allele-defining SNPs (ADSs) within Δ bp of the specified intervals will be discarded")
    parser.add_argument("-d", "--distance", metavar="Δ", default=4, help="minimum distance Δ between ASDs and SNP mask intervals; by default Δ=4")
    return parser

if __name__ == "__main__":
    args = get_parser().parse_args()
    main()
