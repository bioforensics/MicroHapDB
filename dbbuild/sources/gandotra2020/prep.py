#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse
import pandas

def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('outfile')
    return parser


def main(args):
    outdata = {
        'Name': list(),
        'Xref': list(),
        'NumVars': list(),
        'Chrom': list(),
        'OffsetsHg37': list(),
        'OffsetsHg38': list(),
        'VarRef': list(),
    }
    data = pandas.read_csv(args.infile, sep='\t').sort_values(by=['Chrom', 'GRCh38'])
    for marker in data.Marker.unique():
        subdata = data[data.Marker == marker]
        numsnps = len(subdata.GRCh37)
        chrom = f"chr{subdata.Chrom.iloc[0]}"
        offsets37 = [o - 1 for o in subdata.GRCh37]
        offsets38 = [o - 1 for o in subdata.GRCh38]
        offstrs37 = ','.join(map(str, offsets37))
        offstrs38 = ','.join(map(str, offsets38))
        rsids = ','.join([rsid for rsid in subdata.RSID if rsid != '.'])
        newname = marker.replace('KK', 'KKCS').replace('NK', 'KKCS')
        outdata['Name'].append(newname)
        outdata['Xref'].append(None)
        outdata['NumVars'].append(numsnps)
        outdata['Chrom'].append(chrom)
        outdata['OffsetsHg37'].append(offstrs37)
        outdata['OffsetsHg38'].append(offstrs38)
        outdata['VarRef'].append(rsids)
    output = pandas.DataFrame(outdata)
    output.to_csv(args.outfile, sep='\t', index=False)


if __name__ == "__main__":
    main(cli().parse_args())
