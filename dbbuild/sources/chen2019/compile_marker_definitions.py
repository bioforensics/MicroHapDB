#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse
from collections import defaultdict
import pandas


def reformat_markers(markerfile):
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
    markers = markers[markers.MarkerName != 'mh11CP-004']  # Included in ALFRED
    for marker, mdata in markers.groupby('MarkerName'):
        name = 'mh' + marker[2:6] + '-' + marker[6:]
        chrom = 'chr' + mdata.Chrom.iloc[0]
        rsids = mdata.SNPID.tolist()
        offsets = [pos - 1 for pos in mdata.Position]
        rsidstr = ','.join(rsids)
        offsetstr = ','.join(map(str, offsets))
        final_defs['Name'].append(name)
        final_defs['Xref'].append(None)
        final_defs['NumVars'].append(len(rsids))
        final_defs['Chrom'].append(chrom)
        final_defs['OffsetsHg37'].append(None)
        final_defs['OffsetsHg38'].append(offsetstr)
        final_defs['VarRef'].append(rsidstr)
    return pandas.DataFrame(final_defs)


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('markers')
    parser.add_argument('out')
    return parser


def main(args):
    markers = reformat_markers(args.markers)
    markers.to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    main(cli().parse_args())
