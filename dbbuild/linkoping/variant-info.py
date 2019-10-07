#!/usr/bin/env python
import argparse
from collections import defaultdict

cli = argparse.ArgumentParser()
cli.add_argument('linkoping')
cli.add_argument('dbSNP')
args = cli.parse_args()

vardata = dict()
with open(args.dbSNP, 'r') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        chrnum, pos, rsid, ref, alt, *values = line.split('\t')
        chrom = 'chr' + chrnum
        allelestr = ref + ',' + alt
        alleles = ','.join(sorted(allelestr.split(',')))
        vardata[rsid] = (chrom, pos, alleles)

print('Label', 'Chrom', 'Position', 'Alleles', 'rsIDs', 'Source', sep='\t')
with open(args.linkoping, 'r') as fh:
    next(fh)
    for line in fh:
        rsid, label = line.strip().split()
        chrom, pos, alleles = vardata[rsid]
        print(label, chrom, pos, alleles, rsid, 'dbSNP151', sep='\t')
