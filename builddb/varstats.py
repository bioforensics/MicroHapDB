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

from lib import DataSource, BaseMarkerDefinition, CompleteMarkerDefinition, Resolver
from pathlib import Path

markers = list()
for markerfile in Path("sources").glob("*/marker.csv"):
    markers.extend(BaseMarkerDefinition.from_csv(markerfile))
rsids = set()
for marker in markers:
    rsids.update(marker.varref)
resolver = Resolver("databases/dbSNP")
resolver.resolve_rsids(rsids)
mismatches = set()
for marker in markers:
    cm = CompleteMarkerDefinition(marker, resolver)
    print(cm)
    if cm.is_match is False:
        mismatches.add(cm.base.name)
print("Mismatches:", *sorted(mismatches))
