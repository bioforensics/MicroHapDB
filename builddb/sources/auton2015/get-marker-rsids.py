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

import microhapdb
import sys

print("Marker", "RSIDs", sep="\t")
for n, row in microhapdb.markers.iterrows():
    rsids = microhapdb.variantmap[microhapdb.variantmap.Marker == row.Name]

    # Some rsIDs are present in recent versions of dbSNP but not in the 1000
    # Genomes Project Phase 3 VCFs. We need to exclude these when estimating
    # microhap frequencies from the 1KGP data.
    absent_from_1kgp = ["rs772115763", "rs1196416099", "rs377732696", "rs78817707"]
    rsids = rsids[~rsids.Variant.isin(absent_from_1kgp)]

    # Some 1KGP rsIDs have been replaced in more recent versions of dbSNP. We
    # need to revert to the rsID used in the 1KGP for estimating microhap
    # frequencies.
    dbSNP_to_1kgp = {
        "rs73151289": "rs74898010",
        "rs4076758": "rs28970291",
        "rs36190610": "rs72629020",
        "rs10987426": "rs113012024",
        "rs71785313": "rs143830837",
        "rs602427": "rs74812635",
        "rs602875": "rs75324027",
    }
    for newid, oldid in dbSNP_to_1kgp.items():
        rsids = rsids.replace(newid, oldid)

    # Exclude markers for which we don't have all the 1KGP rsIDs.
    nvariants = row.Offsets.count(",") + 1
    if len(rsids) != nvariants:
        msg = "marker {mrkr:s} contains {nv:d} variants, but {nr:d} rsIDs found; skipping".format(
            mrkr=row.Name, nv=nvariants, nr=len(rsids)
        )
        print("[WARNING]", msg, file=sys.stderr)
        continue

    # Yay, print out the marker and all of its component variants.
    print(row.Name, ",".join(rsids.Variant), sep="\t")
