#!/usr/bin/env python

import re
import sys

def parse_block(instream):
    text = instream.read()
    for block in text.split('\f')[1:]:
        block = block.strip()
        block = block.replace(' -', '-')
        if block:
            yield block


def parse_marker(instream):
    for block in parse_block(instream):
        blocklines = re.split(r'\s+', block)
        label = blocklines[0]
        coords = blocklines[1]
        extent = blocklines[2]
        rsids = ','.join(blocklines[3:-2]).replace('nors', 'rs772115763').replace('rs74898010', 'rs73151289').replace('rs28970291', 'rs4076758').replace('rs72629020', 'rs36190610')
        gd = blocklines[-1]
        yield label, coords, extent, rsids, gd


print('Label', 'GRCh37Position', 'Span', 'RSIDs', 'AvgGD', sep='\t')
for data in parse_marker(sys.stdin):
    print(*data, sep='\t')
