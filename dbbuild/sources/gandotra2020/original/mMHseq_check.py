#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

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
