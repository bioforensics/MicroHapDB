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


def collapse_marker_definitions(markers):
    data = {
        'Marker': list(),
        'Offsets37': list(),
        'Offsets38': list(),
        'RSIDs': list(),
    }
    for marker, mdata in markers.groupby('Marker'):
        o37 = [int(v.split(':')[1]) - 1 for v in mdata.Position37 if not pandas.isna(v)]
        o38 = [int(v.split(':')[1]) - 1 for v in mdata.Position38]
        rs = mdata[mdata.RSID != '-'].RSID.tolist()
        data['Marker'].append(marker)
        data['Offsets37'].append(o37)
        data['Offsets38'].append(o38)
        data['RSIDs'].append(rs)
    return pandas.DataFrame(data)


def parse_supp_data(stream):
    for line in stream:
        if not line.startswith('mh'):
            continue
        name, chrom, amplstart, amplend, localoffsets = line.strip().split(',')
        offsets = [int(amplstart) + int(o) for o in localoffsets.split(':')]
        offsets = sorted(offsets)
        if amplstart > amplend:
            amplstart, amplend = amplend, amplstart
        yield name, chrom, int(amplstart), int(amplend), offsets


def check_definitions(markerfile, suppfile):
    finaldefs = {
        'Name': list(),
        'Xref': list(),
        'NumVars': list(),
        'Chrom': list(),
        'OffsetsHg37': list(),
        'OffsetsHg38': list(),
        'VarRef': list(),
    }
    marker_rsids = pandas.read_csv(markerfile, sep='\t')
    markers = collapse_marker_definitions(marker_rsids)
    with open(suppfile, 'r') as fh:
        for name, chrom, astart, aend, offsets in parse_supp_data(fh):
            mdata = markers[markers.Marker == name].iloc[0]
            assert offsets == mdata.Offsets38, mdata
            o37 = ','.join(map(str, mdata.Offsets37))
            if o37 == '':
                o37 = None
            o38 = ','.join(map(str, mdata.Offsets38))
            rsids = ','.join(mdata.RSIDs)
            finaldefs['Name'].append(name)
            finaldefs['Xref'].append(None)
            finaldefs['NumVars'].append(len(offsets))
            finaldefs['Chrom'].append(chrom)
            finaldefs['OffsetsHg37'].append(o37)
            finaldefs['OffsetsHg38'].append(o38)
            finaldefs['VarRef'].append(rsids)
    return pandas.DataFrame(finaldefs)


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('markerrsids')
    parser.add_argument('markerdefs')
    parser.add_argument('out')
    return parser


def main(args):
    markers = check_definitions(args.markerrsids, args.markerdefs)
    markers.to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    main(cli().parse_args())
