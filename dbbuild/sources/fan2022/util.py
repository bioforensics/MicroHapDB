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


Marker = namedtuple("Marker", "name nvars chrom positions varref")


def compile_marker_definitions(infile, outfile):
    markers = list()
    with open(infile, "r") as fh:
        for marker in parse_markers(fh):
            positions = ";".join(map(str, marker.positions))
            rsids = ";".join(marker.varref)
            chrom = f"chr{marker.chrom}"
            entry = (marker.name, None, marker.nvars, "GRCh38", chrom, positions, rsids)
            markers.append(entry)
    colnames = ["Name", "Xref", "NumVars", "Refr", "Chrom", "Positions", "VarRef"]
    table = pd.DataFrame(markers, columns=colnames)
    table.to_csv(outfile, index=False)


def collect_indel_metadata(infile, outfile):
    indels = list()
    with open(infile, "r") as fh:
        for marker, index, refr, alt in parse_indels(fh):
            entry = (marker, index, refr, alt)
            indels.append(entry)
    table = pd.DataFrame(indels, columns=["Marker", "VariantIndex", "Refr", "Alt"])
    table.to_csv(outfile, index=False)


def parse_markers(instream):
    next(instream)
    marker = None
    for line in instream:
        if line.strip() == "":
            continue
        fields = line.rstrip().split("\t")
        pos38 = int(fields[4])
        varref = fields[8] if len(fields) >= 9 else ""
        if fields[0].startswith("MH"):
            if marker is not None:
                yield marker
            name = "mh" + fields[0][2:]
            nvars = int(fields[1])
            chrom = int(fields[2][3:])
            varreflist = list() if varref == "" else [varref]
            marker = Marker(name, nvars, chrom, [pos38], varreflist)
        else:
            marker.positions.append(pos38)
            if varref != "":
                marker.varref.append(varref)
    yield marker


def parse_indels(instream):
    next(instream)
    marker = None
    index = -1
    for line in instream:
        if line.strip() == "":
            continue
        fields = line.rstrip().split("\t")
        if fields[0].startswith("MH"):
            marker = "mh" + fields[0][2:]
            index = -1
        index += 1
        refr = fields[6]
        alts = fields[7].split("/")
        altlengths = sorted([len(alt) for alt in alts])
        if len(refr) > 1 or altlengths[-1] > 1:
            yield marker, index, refr, ";".join(alts)
