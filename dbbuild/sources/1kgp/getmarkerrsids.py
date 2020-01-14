#!/usr/bin/env python
import microhapdb

print('Marker', 'Offsets', 'rsIDs', sep='\t')
for n, row in microhapdb.markers.iterrows():
    result = microhapdb.variantmap[
        (microhapdb.variantmap.Marker == row.Name)
        & (microhapdb.variantmap.Variant.str.startswith('rs'))
    ]
    rsids = list(result.Variant)
    print(row.Name, row.Offsets, ','.join(rsids), sep='\t')
