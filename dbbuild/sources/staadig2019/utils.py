# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import builtins
from collections import defaultdict
from gzip import open as gzopen
from re import findall, search
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


def rename_marker(oldname, chrom):
    return 'mh{chrnum:02d}AT-{serial:s}'.format(chrnum=int(chrom[3:]), serial=oldname[2:4])


def linköping_marker_coords(vcf, mapping):
    rsidcoords = dict()
    for line in vcf:
        if line.startswith('#'):
            continue
        chrnum, posstr, rsid, *values = line.strip().split()
        chrom = chrnum if chrnum.startswith('chr') else 'chr' + chrnum
        pos = int(posstr) - 1
        rsidcoords[rsid] = (chrom, pos)

    marker_rsids = defaultdict(list)
    next(mapping)
    for line in mapping:
        marker, rsid, *values = line.strip().split('\t')
        marker_rsids[marker].append(rsid)

    for marker, rsids in marker_rsids.items():
        positions = [rsidcoords[r][1] for r in rsids]
        assert positions == sorted(positions)
        chrom = rsidcoords[rsids[0]][0]
        name = rename_marker(marker, chrom)
        yield name, chrom, positions, rsids


def linköping_allele_frequencies(freqstream, namestream):
    newnames = dict()
    for line in namestream:
        if line.startswith('Name'):
            continue
        newname, *values = line.split('\t')
        oldname = 'MH' + newname[-2:]
        newnames[oldname] = newname

    for line in freqstream:
        if line.startswith('Microhaplotype'):
            continue
        oldname, allelestr, freqstr = line.strip().split('\t')
        if oldname.endswith(('A', 'B')):
            continue
        name = newnames[oldname]
        allele = ','.join(allelestr)
        yield name, allele, float(freqstr)
