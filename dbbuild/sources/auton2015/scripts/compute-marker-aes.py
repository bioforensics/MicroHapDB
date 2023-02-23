#!/usr/bin/env python
# -------------------------------------------------------------------------------------------------
# Copyright (c) 2019, DHS.
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
import pandas as pd
import sys


def compute_marker_aes(freqfile, popsfile):
    freqs = pd.read_csv(freqfile)
    sample_pops = pd.read_csv(popsfile, sep="\t")
    pops1kgp = set(sample_pops.Population)
    aedata = list()
    for marker, mdata in freqs[freqs.Population.isin(pops1kgp)].groupby("Marker"):
        for population, pdata in mdata.groupby("Population"):
            ae = 1.0 / sum([f**2 for f in pdata.Frequency])
            entry = (marker, population, ae)
            aedata.append(entry)
    table = pd.DataFrame(aedata, columns=["Marker", "Population", "Ae"])
    return table.sort_values(["Marker", "Population"])


def main(args):
    aedata = compute_marker_aes(args.freqs, args.pops)
    aedata.to_csv(sys.stdout, index=False, float_format="%.4f")


def cli():
    parser = ArgumentParser()
    parser.add_argument("freqs")
    parser.add_argument("pops")
    return parser


if __name__ == "__main__":
    main(cli().parse_args())
