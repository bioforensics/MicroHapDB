#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse
import pandas
import rsidx
import sqlite3
import sys


def load_table(filename):
    table = pandas.read_csv(filename, sep="\t")
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
    table[colname] = pandas.Series(newcolumn).astype(int)


def cli():
    desc = (
        "Retrieve SNP positions from dbSNP VCFs using rsIDs, and create a BED file to liftover "
        "SNPs with no RSID"
    )
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("input")
    parser.add_argument("vcf37")
    parser.add_argument("rsidx37")
    parser.add_argument("vcf38")
    parser.add_argument("rsidx38")
    parser.add_argument("output")
    parser.add_argument("novel")
    return parser
    args = cli.parse_args()


def main(args):
    table = load_table(args.input)
    add_coord_column(table, "VCF37", args.rsidx37, args.vcf37)
    add_coord_column(table, "VCF38", args.rsidx38, args.vcf38)
    with open(args.novel, "w") as fh:
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
    table.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main(cli().parse_args())
