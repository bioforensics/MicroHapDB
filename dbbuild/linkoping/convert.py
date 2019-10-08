#!/usr/bin/env python
import sys

print('Variant', 'Microhaplotype', sep='\t')
for line in sys.stdin:
    if not line.startswith('MH'):
        continue
    name, rsid, chrom, pos, refr = line.split('\t')
    label = 'mh{chrnum:02d}AT-{serial:s}'.format(chrnum=int(chrom[3:]), serial=name[2:4])
    print(rsid, label, sep='\t')
