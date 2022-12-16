#!/usr/bin/env python
# -------------------------------------------------------------------------------------------------
# Copyright (c) 2022, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from lib import DataSource
from pathlib import Path

sources = list()
for sourcepath in Path("sources").iterdir():
    if not sourcepath.is_dir():
        continue
    sources.append(DataSource(sourcepath, dbsnp_path="databases/dbSNP"))
sources.sort(key=lambda s: (s.year, s.name.lower()))
mismatches = set()
for source in sources:
    if source.markers is None:
        continue
    for marker in source.markers:
        match = None
        if marker.positions is not None and marker.positions37 is not None:
            match = sorted(marker.positions) == sorted(marker.positions37)
        elif marker.positions is not None and marker.positions38 is not None:
            match = sorted(marker.positions) == sorted(marker.positions38)
        if match is False:
            mismatches.add(marker.name)
        print(marker.name, match, marker.region, marker.region37, marker.region38, sep="\t")
print("Mismatches:", *sorted(mismatches))
