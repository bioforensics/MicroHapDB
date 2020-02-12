#!/usr/bin/env python
import microhapdb
import sys

print('Marker', 'RSIDs', sep='\t')
for n, row in microhapdb.markers.iterrows():
    nvariants = row.Offsets.count(',') + 1
    rsids = microhapdb.variantmap[microhapdb.variantmap.Marker == row.Name]
    absent_from_1kgp = ['rs772115763', 'rs1196416099', 'rs377732696', 'rs78817707']
    rsids = rsids[~rsids.Variant.isin(absent_from_1kgp)]
    dbSNP_to_1kgp = {
        'rs73151289': 'rs74898010',
        'rs4076758': 'rs28970291',
        'rs36190610': 'rs72629020',
        'rs10987426': 'rs113012024',
        'rs71785313': 'rs143830837',
    }
    for newid, oldid in dbSNP_to_1kgp.items():
        rsids = rsids.replace(newid, oldid)
    if len(rsids) != nvariants:
        msg = 'WARNING: marker {:s} contains {:d} variants, but {:d} rsIDs found; skipping'.format(
            row.Name, nvariants, len(rsids)
        )
        print(msg, file=sys.stderr)
        continue
    print(row.Name, ','.join(rsids.Variant), sep='\t')
