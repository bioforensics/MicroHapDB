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

import pandas as pd


def reformat_markers(infile, outtable, outmap):
    definitions = list()
    idmap = list()
    with open(infile, "r") as fh:
        for marker, name, xref, chrom, positions, rsids in parse_marker_table(fh):
            idmap.append((name, marker))
            numvars = len(positions)
            positions = ";".join(map(str, positions))
            rsids = ";".join(rsids)
            entry = (marker, xref, numvars, "GRCh37", chrom, positions, rsids)
            definitions.append(entry)
    colnames = ["Name", "Xref", "NumVars", "Refr", "Chrom", "Positions", "VarRef"]
    deftable = pd.DataFrame(definitions, columns=colnames)
    deftable.to_csv(outtable, index=False)
    maptable = pd.DataFrame(idmap, columns=["InternalName", "Marker"])
    maptable.to_csv(outmap, index=False)


def parse_marker_table(instream):
    marker, name, xref, chrom = [None] * 4
    positions, rsids = list(), list()
    next(instream)
    for line in instream:
        values = line.strip().split()
        if len(values) == 0:
            continue
        elif len(values) == 6:
            if marker is not None:
                yield marker, name, xref, chrom, positions, rsids
            name = values[0]
            positions = list()
            rsids = list()
            marker = values[4]
            xref = values[5]
            chrom = values[2]
            positions.append(int(values[3]))
            rsids.append(values[1])
        else:
            assert len(values) == 3, values
            positions.append(int(values[2]))
            rsids.append(values[0])
    yield marker, name, xref, chrom, positions, rsids


def reformat_frequencies(infile, mapfile, outfile):
    idmaptable = pd.read_csv(mapfile)
    idmap = dict(zip(idmaptable.InternalName, idmaptable.Marker))
    freqs = pd.read_csv(infile, sep="\t").replace(idmap)
    counts_by_marker = dict()
    for markerid, subset in freqs.groupby("MarkerName"):
        counts_by_marker[markerid] = subset.NumberOfObservations.sum()
    freqs.drop(columns=["NumberOfObservations"], inplace=True)
    freqs["Haplotype"] = freqs["Haplotype"].apply(lambda x: "|".join(list(x)))
    freqs["Population"] = "MHDBP-7c055e7ee8"
    freqs.rename(columns={"MarkerName": "Marker", "Haplotype": "Allele"}, inplace=True)
    freqs = freqs[["Marker", "Population", "Allele", "Frequency"]]
    freqs["Count"] = freqs.Marker.apply(lambda x: counts_by_marker[x])
    freqs.to_csv(outfile, index=False, float_format="%.3f")
