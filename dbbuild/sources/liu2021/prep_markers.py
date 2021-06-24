#!/usr/bin/env python
from collections import defaultdict
import sys
marker = None
rsids = defaultdict(list)
for line in sys.stdin:
    fields = line.rstrip().split("\t")
    if fields[0] != "":
        marker = fields[0]
    rsids[marker].append(fields[2])
print("Name", "Xref", "NumVars", "Chrom", "OffsetsHg37", "OffsetsHg38", "VarRef", sep="\t")
for marker, rsidlist in sorted(rsids.items()):
    chrom = marker[2:4]
    if chrom[0] == "0":
        chrom = chrom[1:]
    print(marker, "", len(rsidlist), chrom, "", "", ",".join(rsidlist), sep="\t")
