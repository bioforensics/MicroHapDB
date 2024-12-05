#!/usr/bin/env python

import microhapdb
import pandas as pd
from util import load_markers

markers = load_markers("data/input/marker-0.12.csv", "data/input/marker-aes-0.12.csv", objectify=False)

failed = pd.read_csv("data/intermediate/markers-failed-filter.tsv", sep="\t")
failed_markers = markers[markers.Name.isin(failed.Marker)]
for i, row in failed_markers.iterrows():
    print(row.Chrom, row.Start, row.End, 3, 1, "#990000", 2, sep="\t")

passed = pd.read_csv("data/results/markers-passed-filter.csv")
for i, row in passed.iterrows():
    print(row.Chrom, row.Start, row.End, 3, 1, "#009900", 1, sep="\t")

final = pd.read_csv("data/results/final-panel.tsv", sep="\t")
final_panel = markers[(markers.Name.isin(final.Marker)) | (markers.Chrom == "chrX")]
max_ae = final_panel.Ae.max()
for i, row in final_panel.iterrows():
    shape = 0 if row.Extent < 100 else 1
    width = 0.8 * row.Ae / max_ae
    print(row.Chrom, row.Start, row.End, 1, width, "#0000ff", 1, sep="\t")
