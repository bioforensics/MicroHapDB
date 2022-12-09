#!/usr/bin/env python
#
# -------------------------------------------------------------------------------------------------
# Copyright (c) 2020, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

import glob
import os
import sys

try:
    csvdir = sys.argv[1]
except Exception:
    print("Please specify CSV directory")
    raise SystemExit

path = os.path.join(csvdir, "*.csv")
csvs = glob.glob(path)
markers = set()
populations = set()
for filename in csvs:
    marker, pop, *_ = os.path.basename(filename).split("_")
    markers.add(marker)
    populations.add(pop)
missing = set()
invalid = set()
for marker in sorted(markers):
    for pop in sorted(populations):
        filename = os.path.join(csvdir, f"{marker}_{pop}_freq_table.csv")
        if not os.path.isfile(filename):
            missing.add(filename)
            continue
        with open(filename, "r") as fh:
            data = fh.read().strip()
            if '"Total"' not in data.split("\n")[-1]:
                invalid.add(filename)
print("[Missing files]", *sorted(missing), sep="\n", end="\n\n")
print("[Invalid files]", *sorted(invalid), sep="\n")
