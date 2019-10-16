# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


import builtins
from collections import defaultdict
from gzip import open as gzopen
import sys


def smartopen(filename, mode):
    """Smart file handler

    Determines whether to create a compressed or un-compressed file handle
    based on filename extension.
    """
    if mode not in ('r', 'w'):
        raise ValueError('invalid mode "{}"'.format(mode))
    if filename in ['-', None]:
        filehandle = sys.stdin if mode == 'r' else sys.stdout
        return filehandle
    openfunc = builtins.open
    if filename.endswith('.gz'):
        openfunc = gzopen
        mode += 't'
    return openfunc(filename, mode)


def parse_markers(instream):
    for line in instream:
        if not line.startswith('mh'):
            continue
        name, chrom, amplstart, amplend, localoffsets = line.strip().split(',')
        offsets = [int(amplstart) + int(o) for o in localoffsets.split(':')]
        offsets = sorted(offsets)
        if amplstart > amplend:
            amplstart, amplend = amplend, amplstart
        yield name, chrom, int(amplstart), int(amplend), offsets


def marker_variants(markerstream, vcfstream):
    rsids_by_pos = defaultdict(dict)
    for line in vcfstream:
        if line.startswith('#'):
            continue
        chrom, posstr, rsid, *values = line.split('\t')
        if not chrom.startswith('chr'):
            chrom = 'chr' + chrom
        pos = int(posstr) - 1
        rsids_by_pos[chrom][pos] = rsid

    for name, chrom, astart, aend, offsets in parse_markers(markerstream):
        rsids = list()
        for o in offsets:
            if o in rsids_by_pos[chrom]:
                rsid = rsids_by_pos[chrom][o]
                rsids.append(rsid)
            else:
                # Only variant not present in dbSNP
                assert chrom == 'chr22'
                assert o in (44857884, 44857930, 44857946)
        yield name, chrom, offsets, rsids


def parse_frequencies(markerstream, freqstream):
    alleles = defaultdict(dict)

    markerid = None
    for line in markerstream:
        values = line.strip().split(',')
        if line.startswith('mh'):
            markerid = values[0]
        else:
            haplotype = values[0][2:]
            allelelist = [a.strip() for a in values[1].split(':')]
            allelestr = ','.join(allelelist)
            alleles[markerid][haplotype] = allelestr

    next(freqstream)
    markerid = None
    for line in freqstream:
        values = line.strip().split('\t')
        if line.startswith('mh'):
            markerid = values[0]
        else:
            if markerid == 'mh17PK-86511':
                continue
            haplotype = values[0]
            freqs = values[1:]
            for freq, pop in zip(freqs, ('NL', 'Asia', 'Africa')):
                allele = alleles[markerid][haplotype]
                yield markerid, pop, allele, '{:.4f}'.format(float(freq))
