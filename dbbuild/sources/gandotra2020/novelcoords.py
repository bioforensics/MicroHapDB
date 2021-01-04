#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse
from collections import defaultdict
import pandas


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("bed37")
    parser.add_argument("bed38")
    parser.add_argument("output")
    return parser


def load_mapping(bed37, bed38):
    mapping = defaultdict(dict)
    with open(bed37, "r") as fh37, open(bed38, "r") as fh38:
        for line37, line38 in zip(fh37, fh38):
            chromstr37, _, pos37 = line37.split()
            chromstr38, _, pos38 = line38.split()
            chrom37, chrom38 = chromstr37[3:], chromstr38[3:]
            assert chrom38 == chrom37
            mapping[int(chrom37)][int(pos37)] = int(pos38)
    return mapping


def main(args):
    mapping = load_mapping(args.bed37, args.bed38)
    table = pandas.read_csv(args.input, sep="\t")
    for n, row in table.iterrows():
        if row.RSID == ".":
            table.loc[n, "GRCh38"] = mapping[row.Chrom][row.GRCh37]
    table.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main(cli().parse_args())
