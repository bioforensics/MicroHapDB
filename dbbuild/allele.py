# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from collections import defaultdict, namedtuple
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


def allele_frequencies_alfred(freqfile, indelfile):
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


def allele_frequencies_lovd(locusfile, freqfile):
    alleles = defaultdict(dict)
    with open(locusfile, 'r') as instream:
        locusid = None
        for line in instream:
            values = line.strip().split(',')
            if line.startswith('mh'):
                locusid = values[0]
            else:
                haplotype = values[0][2:]
                allelelist = [a.strip() for a in values[1].split(':')]
                allelestr = ','.join(allelelist)
                alleles[locusid][haplotype] = allelestr

    with open(freqfile, 'r') as instream:
        next(instream)
        locusid = None
        for line in instream:
            values = line.strip().split('\t')
            if line.startswith('mh'):
                locusid = values[0]
            else:
                if locusid == 'mh17PK-86511':
                    continue
                haplotype = values[0]
                freqs = values[1:]
                for freq, pop in zip(freqs, ('NL', 'Asia', 'Africa')):
                    allele = alleles[locusid][haplotype]
                    yield AlleleFreq(locusid, pop, allele, freq)
