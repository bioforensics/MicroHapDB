#!/usr/bin/env python
import microhapdb
import sys

print('Marker', 'RSIDs', sep='\t')
for n, row in microhapdb.markers.iterrows():
    nvariants = row.Offsets.count(',') + 1
    rsids = microhapdb.variantmap[microhapdb.variantmap.Marker == row.Name]
    if len(rsids) != nvariants:
        msg = 'WARNING: marker {:s} contains {:d} variants, but {:d} rsIDs found; skipping'.format(
            row.Name, nvariants, len(rsids)
        )
        print(msg, file=sys.stderr)
    print(row.Name, ','.join(rsids.Variant), sep='\t')
