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
from re import search
from utils import smartopen


def parse_pop_id_mapping(mappingfile):
    m = dict()
    mapping = pandas.read_csv(mappingfile, sep='\t')
    for n, row in mapping.iterrows():
        m[row.ALFRED] = row.ID1KGP
    return m


def parse_freqs(stream, mapping):
    indels = {
        'SI664579L': {'D': 'T', 'I': 'TA'},
        'SI664597L': {'D': 'T', 'I': 'TG'},
        'SI664640A': {'D': 'A', 'I': 'AATAATT'},
    }

    def cleanup(allelestr, markerid, indels):
        allelestr = allelestr.replace('-', ',')
        if 'D' in allelestr or 'I' in allelestr:
            allelestr = allelestr.replace('D', indels[markerid]['D'])
            allelestr = allelestr.replace('I', indels[markerid]['I'])
        return allelestr

    line = next(stream)
    if line.startswith('Created on'):
        next(stream)
    for chunk in stream.read().split('-----------------\n'):
        lines = iter(chunk.split('\n'))
        header1 = next(lines)
        xref, _, _, markerid = header1.split(' | ')
        if '-' not in markerid:
            markerid = markerid[:6] + '-' + markerid[6:]
        header2 = next(lines)
        alleles = header2.split()[2:]
        alleles = [cleanup(a, xref, indels) for a in alleles]
        for line in lines:
            if line.strip() == '':
                continue
            values = line.split('\t')
            popid = search(r'^([^\(]+)\((\S+)\)', values[0]).group(2)
            freqs = values[2:]
            if len(alleles) != len(freqs):
                message = 'WARNING: allele/freq mismatch '
                message += 'for locus ' + markerid
                message += ' and population ' + popid
                message += f'; {len(freqs)} frequencies vs {len(alleles)} alleles'
                raise ValueError(message)
            for allele, freq in zip(alleles, freqs):
                if popid in mapping:
                    popid = mapping[popid]
                yield markerid, popid, allele, float(freq)


def alfred_frequencies(freqstream, mapping):
    freqdata = {
        'Marker': list(),
        'Population': list(),
        'Allele': list(),
        'Frequency': list(),
    }
    for marker, pop, allele, freq in parse_freqs(freqstream, mapping):
        if marker == 'mh05KK-058':
            marker = 'mh15KK-058'
        freqdata['Marker'].append(marker)
        freqdata['Population'].append(pop)
        freqdata['Allele'].append(allele)
        freqdata['Frequency'].append(freq)
    freqs = pandas.DataFrame(freqdata)
    freqs_1kgp = freqs[freqs.Population.isin(mapping.values())]
    freqs = freqs[~freqs.Population.isin(mapping.values())]
    return freqs, freqs_1kgp


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('freqtable')
    parser.add_argument('idmapping')
    parser.add_argument('out')
    parser.add_argument('out1kgp')
    return parser


def main(args):
    mapping = parse_pop_id_mapping(args.idmapping)
    with smartopen(args.freqtable, 'r') as fh:
        freqs, freqs_1kgp = alfred_frequencies(fh, mapping)
    freqs.to_csv(args.out, sep='\t', index=False, float_format='%.04f')
    freqs_1kgp.to_csv(args.out1kgp, sep='\t', index=False, float_format='%.04f')


if __name__ == '__main__':
    main(cli().parse_args())
