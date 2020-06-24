# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


import builtins
from collections import defaultdict
from gzip import open as gzopen
import pandas
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


def marker_variants(markerstream, rsids):
    data = {
        'Marker': list(),
        'Offsets37': list(),
        'Offsets38': list(),
        'RSIDs': list(),
    }
    for marker, mdata in rsids.groupby('Marker'):
        o37 = [int(v.split(':')[1]) - 1 for v in mdata.Position37 if not pandas.isna(v)]
        o38 = [int(v.split(':')[1]) - 1 for v in mdata.Position38]
        rs = mdata[mdata.RSID != '-'].RSID.tolist()
        data['Marker'].append(marker)
        data['Offsets37'].append(o37)
        data['Offsets38'].append(o38)
        data['RSIDs'].append(rs)
    data = pandas.DataFrame(data)

    for name, chrom, astart, aend, offsets in parse_markers(markerstream):
        mdata = data[data.Marker == name].iloc[0]
        o37 = mdata.Offsets37
        o38 = mdata.Offsets38
        assert offsets == o38, mdata
        yield name, len(offsets), chrom, o37, o38, mdata.RSIDs


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
    popids = {
        'Africa': 'MHDBP-3dab7bdd14',
        'Asia': 'MHDBP-936bc36f79',
        'NL': 'MHDBP-383d86606a',
    }
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
                yield markerid, popids[pop], allele, '{:.4f}'.format(float(freq))
