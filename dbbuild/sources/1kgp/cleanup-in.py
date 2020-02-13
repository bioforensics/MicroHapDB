#!/usr/bin/env python
import re
import sys

print('Marker', 'Informativeness', sep='\t')
next(sys.stdin)
for line in sys.stdin:
    if line.startswith('Command:'):
        break
    values = re.split(r'\s+', line.strip())
    marker, inf = values[:2]
    print(marker, '{:.04f}'.format(float(inf)), sep='\t')
