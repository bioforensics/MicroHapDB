# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
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


def marker_coords(vcf, mapping):
    rsidcoords = dict()
    for line in vcf:
        if line.startswith('#'):
            continue
        chrnum, posstr, rsid, *values = line.strip().split()
        chrom = chrnum if chrnum.startswith('chr') else 'chr' + chrnum
        pos = int(posstr) - 1
        rsidcoords[rsid] = (chrom, pos)

    marker_rsids = dict()
    next(mapping)
    for line in mapping:
        marker, pos, span, rsidlist, avggd = line.strip().split('\t')
        marker_rsids[marker] = rsidlist.split(',')

    for marker, rsids in marker_rsids.items():
        positions = [rsidcoords[r][1] for r in rsids]
        rsids = [rsid for pos, rsid in sorted(zip(positions, rsids))]
        positions = sorted(positions)
        chrom = rsidcoords[rsids[0]][0]
        chromlabel = '0X' if chrom == 'chrX' else '{:02d}'.format(int(chrom[3:]))
        mhname = 'mh{chr:s}USC-{n:s}'.format(chr=chromlabel, n=marker)
        yield mhname, chrom, positions, rsids
