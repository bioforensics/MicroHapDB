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
from utils import smartopen


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


def alfred_marker_coords(markerfile, rsidcoords):
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
    for n, row in markers.iterrows():
        name = row.Name
        if name == 'mh05KK-058':
            name = 'mh15KK-058'
        rsids = row.rsIDs.split(',')
        chrom = rsidcoords[rsids[0]][0]
        offsets = [rsidcoords[rsid][1] for rsid in rsids]
        offsetstr = ','.join(map(str, offsets))
        final_defs['Name'].append(name)
        final_defs['Xref'].append(row.Xref)
        final_defs['NumVars'].append(len(offsets))
        final_defs['Chrom'].append(chrom)
        final_defs['OffsetsHg37'].append(None)
        final_defs['OffsetsHg38'].append(offsetstr)
        final_defs['VarRef'].append(row.rsIDs)
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
    markers = alfred_marker_coords(args.markerrsids, rsidcoords)
    markers.to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    main(cli().parse_args())
