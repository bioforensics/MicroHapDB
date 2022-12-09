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

import argparse
from gzip import open as gzopen
import pandas as pd
import subprocess


def compile_marker_definitions(rsids, vcf, outcsv):
    with gzopen(vcf, "rt") as fh:
        rsidcoords = parse_rsid_coords(fh)
    markers = marker_coords(rsids, rsidcoords)
    markers.to_csv(outcsv, index=False)


def parse_rsid_coords(vcf):
    rsidcoords = dict()
    for line in vcf:
        if line.startswith("#"):
            continue
        chrnum, posstr, rsid, *values = line.strip().split()
        chrom = chrnum if chrnum.startswith("chr") else "chr" + chrnum
        pos = int(posstr)
        rsidcoords[rsid] = (chrom, pos)
    return rsidcoords


def marker_coords(markerfile, rsidcoords):
    final_defs = list()
    markers = pd.read_csv(markerfile, sep="\t")
    for n, row in markers.iterrows():
        rsids = row.RSIDs.split(",")
        positions = [rsidcoords[r][1] for r in rsids]
        rsids = [rsid for pos, rsid in sorted(zip(positions, rsids))]
        rsidstr = ";".join(rsids)
        positions = sorted(positions)
        offsetstr = ";".join(map(str, positions))
        chrom = rsidcoords[rsids[0]][0]
        chromlabel = "0X" if chrom == "chrX" else "{:02d}".format(int(chrom[3:]))
        name = f"mh{chromlabel}USC-{row.Label}"
        entry = (name, None, len(positions), "GRCh38", chrom, offsetstr, rsidstr)
        final_defs.append(entry)
    colnames = ["Name", "Xref", "NumVars", "Refr", "Chrom", "Positions", "VarRef"]
    return pd.DataFrame(final_defs, columns=colnames)


def subset_dbsnp(tsv, invcf, rsidx, outvcf):
    markers = pd.read_csv(tsv, sep="\t")
    rsids = list()
    for n, row in markers.iterrows():
        rsids.extend(row.RSIDs.split(","))
    args = ["rsidx", "search", "--header", "--out", outvcf, invcf, rsidx, *rsids]
    subprocess.run(args)


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
