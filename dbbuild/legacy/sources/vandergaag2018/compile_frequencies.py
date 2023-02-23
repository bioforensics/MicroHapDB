#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


import argparse
from collections import defaultdict
import pandas


def parse_haplotypes(stream):
    haplotypes = defaultdict(dict)
    marker = None
    for line in stream:
        values = line.strip().split(',')
        if line.startswith('mh'):
            marker = values[0]
        else:
            hapid = values[0][2:]
            allelelist = [a.strip() for a in values[1].split(':')]
            haplotype = ','.join(allelelist)
            haplotypes[marker][hapid] = haplotype
    return haplotypes


def parse_frequencies(stream, haplotypes):
    freqdata = {
        'Marker': list(),
        'Population': list(),
        'Allele': list(),
        'Frequency': list(),
    }
    next(stream)
    popids = {
        'Africa': 'MHDBP-3dab7bdd14',
        'Asia': 'MHDBP-936bc36f79',
        'NL': 'MHDBP-383d86606a',
    }
    marker = None
    for line in stream:
        values = line.strip().split('\t')
        if line.startswith('mh'):
            marker = values[0]
        else:
            if marker == 'mh17PK-86511':
                continue
            hapid = values[0]
            freqs = values[1:]
            for freq, pop in zip(freqs, ('NL', 'Asia', 'Africa')):
                freqdata['Marker'].append(marker)
                freqdata['Population'].append(popids[pop])
                freqdata['Allele'].append(haplotypes[marker][hapid])
                freqdata['Frequency'].append(float(freq))
    return pandas.DataFrame(freqdata)


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('markerdefs')
    parser.add_argument('freqs')
    parser.add_argument('out')
    return parser


def main(args):
    with open(args.markerdefs, 'r') as fh:
        haplotypes = parse_haplotypes(fh)
    with open(args.freqs, 'r') as fh:
        freqdata = parse_frequencies(fh, haplotypes)
    freqdata.to_csv(args.out, sep='\t', index=False, float_format='%.03f')


if __name__ == '__main__':
    main(cli().parse_args())
