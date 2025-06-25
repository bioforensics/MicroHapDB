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
from matplotlib import pyplot as plt
import pandas as pd


def main(microhaps, codis_strs, outfile):
    mh = pd.read_csv(microhaps, sep="\t")
    codis = pd.read_csv(codis_strs, sep="\t")

    fig = plt.figure(figsize=(6, 4), dpi=300)
    plt.xlim((0, 400))
    plt.ylim((0, 18))
    plt.yticks([0, 3, 6, 9, 12, 15, 18])
    plt.scatter(mh.Extent, mh.Ae, marker=".", color="#1f77b4", label="Proposed 158-plex MH Panel")
    plt.scatter(codis.Extent, codis.Ae, marker="v", color="#ff7f0e", label="CODIS STR Panel")
    plt.xlabel("Extent (bp)")
    plt.ylabel("Effective Number of Alleles ($A_e$)")
    ax = plt.gca()
    ax.set_axisbelow(True)
    plt.grid(True, linestyle="--", color="#eeeeee")
    ax.legend(loc="upper left")
    fig.savefig(outfile, bbox_inches="tight")


def get_parser():
    parser = ArgumentParser()
    parser.add_argument("microhaps", help="path to final microhap panel in CSV format")
    parser.add_argument("codis", help="path to CODIS STR panel in CSV format")
    parser.add_argument("outfile", help="output filename")
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args.microhaps, args.codis, args.outfile)
