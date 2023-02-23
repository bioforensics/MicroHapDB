#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse
import pandas


def marker_names(markerfile):
    newnames = dict()
    markers = pandas.read_csv(markerfile, sep='\t')
    for n, row in markers.iterrows():
        oldname = 'MH' + row.Name[-2:]
        newnames[oldname] = row.Name
    return newnames


def linköping_allele_frequencies(freqfile, names):
    freqdata = {
        'Marker': list(),
        'Population': list(),
        'Allele': list(),
        'Frequency': list(),
    }
    freqs = pandas.read_csv(freqfile, sep='\t')
    freqs = freqs[~freqs['Microhaplotype name'].str.endswith(('A', 'B'))]
    for n, row in freqs.iterrows():
        name = names[row['Microhaplotype name']]
        haplotype = ','.join(row['Haplotype variant'])
        freqdata['Marker'].append(name)
        freqdata['Population'].append('MHDBP-7c055e7ee8')
        freqdata['Allele'].append(haplotype)
        freqdata['Frequency'].append(row['Frequency'])
    return pandas.DataFrame(freqdata)


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('markerfile')
    parser.add_argument('freqfile')
    parser.add_argument('out')
    return parser


def main(args):
    names = marker_names(args.markerfile)
    freqs = linköping_allele_frequencies(args.freqfile, names)
    freqs.to_csv(args.out, sep='\t', index=False, float_format='%.04f')


if __name__ == '__main__':
    main(cli().parse_args())
