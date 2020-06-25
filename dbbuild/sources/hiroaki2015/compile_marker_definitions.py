#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse
import builtins
from collections import defaultdict
from gzip import open as gzopen
import pandas


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


def parse_rsid_coords(vcf):
    rsidcoords = dict()
    for line in vcf:
        if line.startswith('#'):
            continue
        chrnum, posstr, rsid, *values = line.strip().split()
        chrom = chrnum if chrnum.startswith('chr') else 'chr' + chrnum
        pos = int(posstr) - 1
        rsidcoords[rsid] = (chrom, pos)
    return rsidcoords


def marker_coords(markerfile, rsidcoords):
    final_defs = {
        'Name': list(),
        'Xref': list(),
        'NumVars': list(),
        'Chrom': list(),
        'OffsetsHg37': list(),
        'OffsetsHg38': list(),
        'VarRef': list(),
    }
    markers = pandas.read_csv(markerfile, sep='\t')
    markers = markers[markers.Marker != 'mh04NH-08']  # Problematic marker definition; see README
    for n, row in markers.iterrows():
        rsids = row.Variants.split(',')
        offsets = [rsidcoords[r][1] for r in rsids]
        offsetstr = ','.join(map(str, offsets))
        final_defs['Name'].append(row.Marker)
        final_defs['Xref'].append(None)
        final_defs['NumVars'].append(len(offsets))
        final_defs['Chrom'].append(row.Chrom)
        final_defs['OffsetsHg37'].append(None)
        final_defs['OffsetsHg38'].append(offsetstr)
        final_defs['VarRef'].append(row.Variants)
    return pandas.DataFrame(final_defs)


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('markerrsids')
    parser.add_argument('vcf')
    parser.add_argument('out')
    return parser


def main(args):
    with smartopen(args.vcf, 'r') as fh:
        rsidcoords = parse_rsid_coords(fh)
    markers = marker_coords(args.markerrsids, rsidcoords)
    markers.to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    main(cli().parse_args())
