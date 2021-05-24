#!/usr/bin/env python

import argparse
import pandas
import sys

cli = argparse.ArgumentParser()
cli.add_argument("ids")
cli.add_argument("table")
args = cli.parse_args()

with open(args.ids, "r") as fh:
    idlist = fh.read().splitlines()
data = pandas.read_csv(args.table, sep="\t")
data = data[data.Name.isin(idlist)]
data["OffsetsHg37"] = None
data["OffsetsHg38"] = None
data["NumVars"] = data.VarRef.apply(lambda x: x.count(",") + 1)
data["Name"] = data.Name.apply(lambda x: x.replace("KKCS", "KK") + ".db")
data.to_csv(sys.stdout, sep="\t", index=False)
