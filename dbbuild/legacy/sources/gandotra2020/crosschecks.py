#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse
import pandas
from pyfaidx import Fasta as FastaIdx
import re
import sys


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("fasta37")
    parser.add_argument("fasta38")
    return parser


def compare_marker_sequences(data, fasta37, fasta38):
    index37 = FastaIdx(fasta37)
    index38 = FastaIdx(fasta38)
    markers = sorted(data.Marker.unique())
    mismatches = 0
    for marker in markers:
        subdata = data[(data.Marker == marker) & (data.VCF37 != 0)]
        chrom = f"chr{subdata.Chrom.iloc[0]}"
        min37, max37 = subdata.VCF37.min() - 1, subdata.VCF37.max()
        min38, max38 = subdata.VCF38.min() - 1, subdata.VCF38.max()
        seq37 = str(index37[chrom][min37:max37]).upper()
        seq38 = str(index38[chrom][min38:max38]).upper()
        if seq37 != seq38:
            mismatches += 1
            print(
                "Mismatched sequences:",
                marker,
                f"hg37={chrom37}:{min37}-{max37}",
                f"hg38={chrom38}:{min38}-{max38}",
                f"\n  {seq37}",
                f"\n  {seq38}",
                file=sys.stderr,
            )
    if mismatches > 0:
        raise SystemError()


def compare_snp_offsets(data):
    markers = sorted(data.Marker.unique())
    for marker in markers:
        subdata = data[(data.Marker == marker)]
        min37, min38 = subdata.GRCh37.min(), subdata.GRCh38.min()
        offsets37 = sorted([pos - min37 for pos in subdata.GRCh37])
        offsets38 = sorted([pos - min38 for pos in subdata.GRCh38])
        if offsets37 != offsets38:
            print(f"{marker}\n  {offsets37}\n  {offsets38}")


def main(args):
    data = pandas.read_csv(args.input, sep="\t")
    compare_marker_sequences(data, args.fasta37, args.fasta38)
    compare_snp_offsets(data)


if __name__ == "__main__":
    main(cli().parse_args())
