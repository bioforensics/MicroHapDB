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

import pandas as pd


def compile_marker_definitions(infile, outfile):
    markers = marker_coords(infile)
    markers.to_csv(outfile, index=False)


def marker_coords(infile):
    definitions = list()
    markers = pd.read_csv(infile, sep="\t")
    for n, row in markers.iterrows():
        rsids = row.RSIDs.replace(",", ";")
        numvars = rsids.count(";") + 1
        chrom = row.GRCh37Position.split(":")[0]
        chromlabel = "0X" if chrom == "X" else "{:02d}".format(int(chrom))
        name = f"mh{chromlabel}USC-{row.Label}"
        entry = (name, None, numvars, None, f"chr{chrom}", None, rsids)
        definitions.append(entry)
    colnames = ["Name", "Xref", "NumVars", "Refr", "Chrom", "Positions", "VarRef"]
    return pd.DataFrame(definitions, columns=colnames)


def text2table(infile, outfile):
    with open(infile, "r") as infh, open(outfile, "w") as outfh:
        print("Label", "GRCh37Position", "Span", "RSIDs", "AvgGD", sep="\t", file=outfh)
        for data in parse_marker(infh):
            print(*data, sep="\t", file=outfh)


def parse_marker(instream):
    for block in parse_block(instream):
        blocklines = re.split(r"\s+", block)
        label = blocklines[0]
        coords = blocklines[1]
        extent = blocklines[2]
        rsids = (
            ",".join(blocklines[3:-2])
            .replace("nors", "rs772115763")
            .replace("rs74898010", "rs73151289")
            .replace("rs28970291", "rs4076758")
            .replace("rs72629020", "rs36190610")
        )
        gd = blocklines[-1]
        yield label, coords, extent, rsids, gd


def parse_block(instream):
    text = instream.read()
    for block in text.split("\f")[1:]:
        block = block.strip()
        block = block.replace(" -", "-")
        if block:
            yield block
