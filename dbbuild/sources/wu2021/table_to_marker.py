#!/usr/bin/env python

import pandas as pd

table = pd.read_csv("tables2.txt", sep="\t", skiprows=3, header=None)
newtable = list()
for i, row in table.iterrows():
    newrow = (
        row[1],
        None,
        row[7],
        row[0],
        None,
        None,
        row[8],
    )
    newtable.append(newrow)
colnames = ["Name", "Xref", "NumVars", "Chrom", "OffsetsHg37", "OffsetsHg38", "VarRef"]
markers = pd.DataFrame(newtable, columns=colnames)
markers.to_csv("marker-orig.tsv", sep="\t", index=False)
