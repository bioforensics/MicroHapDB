#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2021, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse
import glob
import os
import pandas
import re
import sys

def parse_filename(filename):
    marker, population, *_ = os.path.basename(filename).split("_")
    marker = re.sub(r"-", "CS-", marker)
    population = f"mMHseq-{population}"
    return marker, population

cli = argparse.ArgumentParser()
cli.add_argument("csvdir")
args = cli.parse_args()

datalist = list()
for csvfile in glob.glob(f"{args.csvdir}/*.csv"):
    marker, population = parse_filename(csvfile)
    data = pandas.read_csv(csvfile, skiprows=2, skipfooter=1, engine="python", header=None)
    data["Allele"] = data[0].apply(lambda x: x.replace("-", ","))
    data["Marker"] = marker
    data["Population"] = population
    data = data.rename(columns={2: "Frequency"})[["Marker", "Population", "Allele", "Frequency"]]
    datalist.append(data)

pandas.concat(datalist).round(4).to_csv(sys.stdout, sep="\t", index=False)
