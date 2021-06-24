#!/usr/bin/env python
import argparse
import pandas
import rsidx
import sqlite3
import sys

cli = argparse.ArgumentParser()
cli.add_argument("tsv")
cli.add_argument("vcf")
cli.add_argument("rsidx")
args = cli.parse_args()

markers = pandas.read_csv(args.tsv, sep="\t")
indel_rsids = [varlist[0] for varlist in markers.VarRef.str.split(",")]
markers_by_rsid = dict()
for i, row in markers.iterrows():
    rsids = row.VarRef.split(",")
    for rsid in rsids:
        markers_by_rsid[rsid] = row.Name

print("Marker", "VariantIndex", "Refr", "Alt", sep="\t")
with sqlite3.connect(args.rsidx) as dbconn:
    for line in rsidx.search.search(indel_rsids, dbconn, args.vcf):
        fields = line.split("\t")
        rsid = fields[2]
        marker = markers_by_rsid[rsid]
        print(marker, 0, fields[3], fields[4], sep="\t")
