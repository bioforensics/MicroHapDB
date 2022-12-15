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

from collections import namedtuple
import pandas as pd


MarkerInfo = namedtuple("MarkerInfo", "name chrom nsnps off37 rsids")


def format_markers(infile, outfile):
    table = list()
    with open(infile, "r") as fh:
        for marker in parse_markers(fh):
            offsets = ";".join(map(str, marker.off37))
            rsids = ";".join(marker.rsids)
            table.append((marker.name, None, marker.nsnps, "GRCh37", marker.chrom, offsets, rsids))
    colnames = ["Name", "Xref", "NumVars", "Refr", "Chrom", "Positions", "VarRef"]
    data = pd.DataFrame(table, columns=colnames)
    data.to_csv(outfile, index=False)


def parse_markers(infile):
    marker = None
    next(infile)
    for line in infile:
        line = line.rstrip()
        if line == "":
            continue
        fields = line.split()
        if len(fields) == 8:
            if marker is not None:
                assert marker.nsnps == len(marker.off37)
                yield marker
            offset = int(fields[4]) - 1
            rsid = fields[5]
            marker = MarkerInfo(
                name=fields[0],
                chrom=f"chr{fields[1]}",
                nsnps=int(fields[3]),
                off37=[offset],
                rsids=[],
            )
            if rsid != ".":
                marker.rsids.append(rsid)
        else:
            offset = int(fields[0]) - 1
            rsid = fields[1]
            marker.off37.append(offset)
            if rsid != ".":
                marker.rsids.append(rsid)
    assert marker.nsnps == len(marker.off37), marker.name
    yield marker
