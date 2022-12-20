#!/usr/bin/env python

import pandas as pd

table = pd.read_csv("tables2.txt", sep="\t", skiprows=3, header=None)
newtable = list()
for i, row in table.iterrows():
    rsids = row[8].replace(",", ";")
    rsids = rsids.replace("rs74812635", "rs602427").replace("rs75324027", "rs602875")
    newrow = (
        row[1],
        None,
        row[7],
        None,
        row[0],
        None,
        rsids,
    )
    newtable.append(newrow)
colnames = ["Name", "Xref", "NumVars", "Refr", "Chrom", "Positions", "VarRef"]
markers = pd.DataFrame(newtable, columns=colnames)
markers.to_csv("marker.csv", index=False)
