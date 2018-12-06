# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from collections import namedtuple
import json
from re import search


AlleleFreq = namedtuple('AlleleFreq', 'locusid, popid, allele, freq')


def load_indel_alleles(infile):
    with open(infile, 'r') as instream:
        return json.load(instream)


def cleanup(allelestr, locusid, indels):
    allelestr = allelestr.replace('-', ',')
    if 'D' in allelestr or 'I' in allelestr:
        allelestr = allelestr.replace('D', indels[locusid]['D'])
        allelestr = allelestr.replace('I', indels[locusid]['I'])
    return allelestr


def allele_frequencies(freqfile, indelfile):
    indels = load_indel_alleles(indelfile)
    with open(freqfile, 'r') as instream:
        line = next(instream)
        if line.startswith('Created on'):
            next(instream)
        chunks = instream.read().split('-----------------\n')

    for chunk in chunks:
        lines = iter(chunk.split('\n'))
        header1 = next(lines)
        locusid = header1.split()[0]
        header2 = next(lines)
        alleles = header2.split()[2:]
        alleles = [cleanup(a, locusid, indels) for a in alleles]
        for line in lines:
            if line.strip() == '':
                continue
            values = line.split('\t')
            popid = search(r'^([^\(]+)\((\S+)\)', values[0]).group(2)
            freqs = values[2:]
            if len(alleles) != len(freqs):
                message = 'WARNING: allele/freq mismatch '
                message += 'for locus ' + locusid
                message += ' and population ' + popid
                message += '; {nfreq} frequencies vs {nall} alleles'.format(
                    nfreq=len(freqs), nall=len(alleles)
                )
                raise ValueError(message)
            for allele, freq in zip(alleles, freqs):
                yield AlleleFreq(locusid, popid, allele, freq)
