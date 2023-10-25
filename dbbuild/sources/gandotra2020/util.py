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

from collections import defaultdict
from glob import glob
import os
import pandas as pd
from pyfaidx import Fasta as FastaIdx
import re
import rsidx
import sqlite3
import sys


def prepare_markers(infile, outfile):
    outdata = list()
    data = pd.read_csv(infile, sep="\t").sort_values(by=["Chrom", "GRCh38"])
    for marker in data.Marker.unique():
        subdata = data[data.Marker == marker]
        numsnps = len(subdata.GRCh38)
        chrom = f"chr{subdata.Chrom.iloc[0]}"
        offstrs38 = ";".join(map(str, subdata.GRCh38))
        rsids = ";".join([rsid for rsid in subdata.RSID if rsid != "."])
        entry = (marker, None, numsnps, "GRCh38", chrom, offstrs38, rsids)
        outdata.append(entry)
    colnames = ["Name", "Xref", "NumVars", "Refr", "Chrom", "Positions", "VarRef"]
    output = pd.DataFrame(outdata, columns=colnames)
    output.to_csv(outfile, index=False)


def compile_freqs(csvdir, outcsv):
    datalist = list()
    for csvfile in glob(f"{csvdir}/*.csv"):
        marker, population = parse_csv_filename(csvfile)
        data = pd.read_csv(csvfile, skiprows=2, skipfooter=1, engine="python", header=None)
        data["Allele"] = data[0].apply(lambda x: x.replace("-", "|"))
        data["Count"] = data[1].sum()
        data["Marker"] = marker
        data["Population"] = population
        columns = ["Marker", "Population", "Allele", "Frequency", "Count"]
        data = data.rename(columns={2: "Frequency"})[columns]
        datalist.append(data)
    pd.concat(datalist).round(4).to_csv(outcsv, index=False, float_format="%.4f")


def parse_csv_filename(filename):
    marker, population, *_ = os.path.basename(filename).split("_")
    population = f"mMHseq-{population}"
    return marker, population


def novel_coords(bed37, bed38, infile, outfile):
    mapping = load_mapping(bed37, bed38)
    table = pd.read_csv(infile, sep="\t")
    for n, row in table.iterrows():
        if row.RSID == ".":
            table.loc[n, "GRCh38"] = mapping[row.Chrom][row.GRCh37]
    table.to_csv(outfile, sep="\t", index=False)


def load_mapping(bed37, bed38):
    mapping = defaultdict(dict)
    with open(bed37, "r") as fh37, open(bed38, "r") as fh38:
        for line37, line38 in zip(fh37, fh38):
            chromstr37, _, pos37 = line37.split()
            chromstr38, _, pos38 = line38.split()
            chrom37, chrom38 = chromstr37[3:], chromstr38[3:]
            assert chrom38 == chrom37
            mapping[int(chrom37)][int(pos37)] = int(pos38)
    return mapping


def cross_checks(infile, fasta37, fasta38):
    data = pd.read_csv(infile, sep="\t")
    compare_marker_sequences(data, fasta37, fasta38)
    compare_snp_offsets(data)


def compare_marker_sequences(data, fasta37, fasta38):
    index37 = FastaIdx(fasta37)
    index38 = FastaIdx(fasta38)
    markers = sorted(data.Marker.unique())
    mismatches = 0
    for marker in markers:
        subdata = data[(data.Marker == marker) & (data.VCF37 != 0)]
        chrom = f"chr{subdata.Chrom.iloc[0]}"
        min37, max37 = subdata.VCF37.min() - 1, subdata.VCF37.max()
        min38, max38 = subdata.VCF38.min() - 1, subdata.VCF38.max()
        seq37 = str(index37[chrom][min37:max37]).upper()
        seq38 = str(index38[chrom][min38:max38]).upper()
        if seq37 != seq38:
            mismatches += 1
            print(
                "Mismatched sequences:",
                marker,
                f"hg37={chrom37}:{min37}-{max37}",
                f"hg38={chrom38}:{min38}-{max38}",
                f"\n  {seq37}",
                f"\n  {seq38}",
                file=sys.stderr,
            )
    if mismatches > 0:
        raise SystemError()


def compare_snp_offsets(data):
    markers = sorted(data.Marker.unique())
    for marker in markers:
        subdata = data[(data.Marker == marker)]
        min37, min38 = subdata.GRCh37.min(), subdata.GRCh38.min()
        offsets37 = sorted([pos - min37 for pos in subdata.GRCh37])
        offsets38 = sorted([pos - min38 for pos in subdata.GRCh38])
        if offsets37 != offsets38:
            print(f"{marker}\n  {offsets37}\n  {offsets38}")


def vcfadd(infile, vcf37, rsidx37, vcf38, rsidx38, outtsv, outbed):
    table = load_table(infile)
    add_coord_column(table, "VCF37", rsidx37, vcf37)
    add_coord_column(table, "VCF38", rsidx38, vcf38)
    with open(outbed, "w") as fh:
        for n, row in table.iterrows():
            if row.RSID == ".":
                print(f"chr{row.Chrom} {row.GRCh37 - 1} {row.GRCh37}", file=fh)
                continue
            if row.GRCh37 != row.VCF37 or row.GRCh38 != row.VCF38:
                print(
                    "Mismatch",
                    row.Marker,
                    row.RSID,
                    row.GRCh37,
                    row.VCF37,
                    row.GRCh37 == row.VCF37,
                    row.GRCh38,
                    row.VCF38,
                    row.GRCh38 == row.VCF38,
                    file=sys.stderr,
                )
    table.to_csv(outtsv, sep="\t", index=False)


def load_table(filename):
    table = pd.read_csv(filename, sep="\t")
    table["GRCh37"] = table.GRCh37.replace(".", 0).astype(int)
    table["GRCh38"] = table.GRCh38.replace(".", 0).astype(int)
    return table


def add_coord_column(table, colname, rsidxdb, vcf):
    positions = dict()
    rsids = list(table[table.RSID != "."].RSID)
    with sqlite3.connect(rsidxdb) as dbconn:
        for line in rsidx.search.search(rsids, dbconn, vcf):
            chrom, position, rsid, *values = line.strip().split()
            positions[rsid] = position
    newcolumn = list()
    for n, row in table.iterrows():
        if row.RSID == ".":
            newcolumn.append(0)
        else:
            position = int(positions[row.RSID])
            newcolumn.append(position)
    table[colname] = pd.Series(newcolumn).astype(int)
