#!/usr/bin/env python
from argparse import ArgumentParser, FileType

cli = ArgumentParser()
cli.add_argument('freqs', type=FileType('r'))
cli.add_argument('defs', type=FileType('r'))
args = cli.parse_args()

labels = dict()
next(args.defs)
for line in args.defs:
    rsid, label = line.strip().split('\t')
    suffix = label[-2:]
    labels[suffix] = label

print('Label', 'Allele', 'Frequency', sep='\t')
next(args.freqs)
for line in args.freqs:
    name, allelestr, freq = line.strip().split('\t')
    if name.endswith(('A', 'B')):
        continue
    suffix = name[2:4]
    label = labels[suffix]
    allele = ','.join(allelestr)
    print(label, allele, freq, sep='\t')
