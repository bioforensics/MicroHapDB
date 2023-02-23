#!/usr/bin/env python
# -------------------------------------------------------------------------------------------------
# Copyright (c) 2019, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

import re
import sys

print("Marker", "Informativeness", sep="\t")
next(sys.stdin)
for line in sys.stdin:
    if line.startswith("Command:"):
        break
    values = re.split(r"\s+", line.strip())
    marker, inf = values[:2]
    print(marker, "{:.04f}".format(float(inf)), sep=",")
