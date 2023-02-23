#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse
from collections import defaultdict
import pandas


def reformat_frequencies(freqfile):
    freqdata = {
        'Marker': list(),
        'Population': list(),
        'Allele': list(),
        'Frequency': list(),
    }
    freqs = pandas.read_csv(freqfile, sep='\t')
    for n, row in freqs.iterrows():
        standardname = 'mh' + row.MarkerName[2:6] + '-' + row.MarkerName[6:]
        haplotype = ','.join(row.Allele)
        freqdata['Marker'].append(standardname)
        freqdata['Population'].append('MHDBP-48c2cfb2aa')
        freqdata['Allele'].append(haplotype)
        freqdata['Frequency'].append(row.Frequency)
    return pandas.DataFrame(freqdata)


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('freqs')
    parser.add_argument('out')
    return parser


def main(args):
    freqs = reformat_frequencies(args.freqs)
    freqs.to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    main(cli().parse_args())
