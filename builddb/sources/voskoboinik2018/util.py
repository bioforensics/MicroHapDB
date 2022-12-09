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
from pathlib import Path
import re
import subprocess
import sys


def compile_marker_definitions(dbsnp38, rsidx38, dir1kgp, table1):
    markernames = list()
    for name, region in parse_marker_definitions(table1):
        markernames.append(name)
        get_variants(name, region, dir1kgp)
        filter_variants(name)
        transfer_variants(name, dbsnp38, rsidx38, refr="GRCh38")
    finalize_marker_definitions(markernames)


def parse_marker_definitions(table1file):
    data = pd.read_csv(table1file, sep="\t")
    data["End"] = data["GRCh37Pos"] + data["Length"]
    for n, row in data.iterrows():
        name = f"mh{row.Chrom:02d}LV-{n+1:02d}"
        region = f"{row.Chrom}:{row.GRCh37Pos}-{row.End}"
        yield name, region


def get_variants(markername, region, path):
    chrom = int(markername[2:4])
    infile = Path(path) / f"chr{chrom}.vcf.gz"
    outfile = f"{markername}-unfiltered.vcf"
    command = ["tabix", infile, region]
    with open(outfile, "w") as fh:
        subprocess.check_call(command, universal_newlines=True, stdout=fh)


def filter_variants(markername):
    infile = f"{markername}-unfiltered.vcf"
    outfile = f"{markername}-filtered.vcf"
    variants = dict()
    with open(infile, "r") as infh, open(outfile, "w") as outfh:
        for var in parse_vcf(infh):
            if not var.is_snp:
                continue
            if var.is_rare():
                continue
            if var.af_slug in variants:
                print(
                    "DUPLICATE", variants[var.af_slug].rsid, var.rsid, var.af_slug, file=sys.stderr
                )
            else:
                variants[var.af_slug] = var
        variants_to_keep = sorted(variants.values(), key=lambda v: v.location)
        rsidmap = {"rs113012024": "rs10987426"}
        for var in variants_to_keep:
            if var.rsid in rsidmap:
                var._fields[2] = rsidmap[var.rsid]
            print(var, file=outfh)


def parse_vcf(infile):
    for line in infile:
        if line.startswith("#"):
            continue
        yield Variant(line)


def transfer_variants(markername, vcf, rsidx, refr="GRCh38"):
    infile = f"{markername}-filtered.vcf"
    rsids = list()
    with open(infile, "r") as fh:
        for line in fh:
            if line.strip() == "" or line.startswith("#"):
                continue
            values = line.split("\t")
            rsid = values[2]
            rsids.append(rsid)
    outfile = f"{markername}-{refr}.vcf"
    command = ["rsidx", "search", vcf, rsidx] + rsids
    with open(outfile, "w") as fh:
        subprocess.check_call(command, universal_newlines=True, stdout=fh)


def finalize_marker_definitions(markernames):
    data = list()
    for name in markernames:
        offsets = list()
        with open(f"{name}-GRCh38.vcf", "r") as fh:
            chrom, rsids38, offsets38 = vcf_to_rsid_offsets(fh)
            offsets = ";".join(map(str, offsets38))
            rsids = ";".join(rsids38)
        entry = (name, None, len(rsids38), "GRCh38", chrom, offsets, rsids)
        data.append(entry)
    colnames = ["Name", "Xref", "NumVars", "Refr", "Chrom", "Positions", "VarRef"]
    df = pd.DataFrame(data, columns=colnames)
    df.to_csv("marker.csv", index=False)


def vcf_to_rsid_offsets(vcfstream):
    offsets = list()
    rsids = list()
    for line in vcfstream:
        chrom, pos, rsid, *values = line.strip().split("\t")
        markerchrom = "chr" + chrom
        offsets.append(int(pos) - 1)
        rsids.append(rsid)
    return markerchrom, rsids, offsets


class Variant:
    def __init__(self, line):
        self._fields = line.strip().split("\t")

    def __str__(self):
        return "\t".join(self._fields)

    @property
    def rsid(self):
        return self._fields[2]

    @property
    def location(self):
        return int(self._fields[0]), int(self._fields[1])

    @property
    def af_slug(self):
        slug = list()
        for kvp in self._fields[7].split(";"):
            if kvp.count("=") != 1:
                continue
            key, value = kvp.split("=")
            if key.endswith("AF"):
                slug.append(key + "=" + value)
        return ";".join(slug)

    @property
    def is_snp(self):
        for allele in self._fields[3].split(",") + self._fields[4].split(","):
            if len(allele) > 1:
                return False
        return True

    def is_rare(self, mincount=2):
        countmatch = re.search(r"AC=(\d+)", self._fields[7])
        assert countmatch
        count = int(countmatch.group(1))
        if count < mincount:
            return True
        return False
