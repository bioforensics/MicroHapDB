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


def rename_marker(oldname, chrom):
    return 'mh{chrnum:02d}AT-{serial:s}'.format(chrnum=int(chrom[3:]), serial=oldname[2:4])


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


def linköping_marker_coords(markerfile, rsidcoords):
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
    markers['Microhaplotype name'] = markers['Microhaplotype name'].apply(lambda x: x[:-1] if x.endswith(('A', 'B')) else x)
    for marker, mdata in markers.groupby('Microhaplotype name'):
        name = rename_marker(marker, mdata.Chromosome.iloc[0])
        rsids = mdata['rs-number'].tolist()
        offsets37 = [pos - 1 for pos in mdata.Position.tolist()]
        offsets38 = [rsidcoords[r][1] for r in rsids]
        rsidstr = ','.join(rsids)
        offsetstr37 = ','.join(map(str, offsets37))
        offsetstr38 = ','.join(map(str, offsets38))
        final_defs['Name'].append(name)
        final_defs['Xref'].append(None)
        final_defs['NumVars'].append(len(offsets38))
        final_defs['Chrom'].append(mdata.Chromosome.iloc[0])
        final_defs['OffsetsHg37'].append(offsetstr37)
        final_defs['OffsetsHg38'].append(offsetstr38)
        final_defs['VarRef'].append(rsidstr)
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
    markers = linköping_marker_coords(args.markerrsids, rsidcoords)
    markers.to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    main(cli().parse_args())
