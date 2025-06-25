#!/usr/bin/env python

# -------------------------------------------------------------------------------------------------
# Copyright (c) 2025, DHS.
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
import matplotlib
from matplotlib import pyplot as plt
from util import load_markers


def main(microhaps, outfile):
    fig = plt.figure(figsize=(6, 4), dpi=300)
    microhaps = bin_by_ads(microhaps)
    plot_by_bins(microhaps)
    plt.xlim((0, 275))
    plt.ylim((0, 18))
    plt.yticks([0, 3, 6, 9, 12, 15, 18])
    plt.xlabel("Extent (bp)")
    plt.ylabel("Effective Number of Alleles ($A_e$)")
    ax = plt.gca()
    ax.set_axisbelow(True)
    plt.grid(True, linestyle="--", color="#eeeeee")
    ax.legend(loc="upper left", fontsize=8)
    fig.savefig(outfile, bbox_inches="tight")


def bin_by_ads(microhaps):
    microhaps = microhaps[microhaps.NumVars > 1].copy()
    microhaps["BinByADS"] = microhaps["NumVars"].apply(lambda nv: nv if nv % 2 == 0 else nv - 1)
    microhaps["BinByADS"] = microhaps["BinByADS"].apply(lambda ads: ads if ads < 10 else 10)
    return microhaps


def plot_by_bins(microhaps):
    set1 = matplotlib.colormaps["Set1"]
    max_num_snps = microhaps.NumVars.max()
    for numsnps, data in microhaps.groupby("BinByADS"):
        colorindex = int((numsnps / 2) - 1)
        color = set1(colorindex)
        markersizes = [nv / max_num_snps * 15 for nv in data.NumVars]
        label = f"{numsnps}-{numsnps+1} ADSs" if numsnps < 10 else "≥10 ADSs"
        marker = "." if numsnps == 2 else "o"
        plt.scatter(data.Extent, data.Ae, marker=marker, s=markersizes, color=color, alpha=0.8, label=label)


def get_parser():
    parser = ArgumentParser()
    parser.add_argument("microhaps", help="path to final microhap panel in CSV format")
    parser.add_argument("aes", help="path to MicroHapDB Ae table in CSV format")
    parser.add_argument("outfile", help="output filename")
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    markers = load_markers(args.microhaps, args.aes, objectify=False)
    main(markers, args.outfile)
