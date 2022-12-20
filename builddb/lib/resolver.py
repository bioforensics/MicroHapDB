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

from pathlib import Path
import rsidx
import sqlite3

class Resolver:
    def __init__(self, dbsnp_path):
        self.dbsnp_path = Path(dbsnp_path)
        self.coords_by_rsid = dict(
            GRCh37=dict(),
            GRCh38=dict(),
        )

    def resolve_rsids(self, rsids):
        vcf = self.dbsnp_path / "dbSNP_GRCh37.vcf.gz"
        idx = self.dbsnp_path / "dbSNP_GRCh37.rsidx"
        Resolver.rsidx_search(rsids, vcf, idx, self.coords_by_rsid["GRCh37"])
        vcf = self.dbsnp_path / "dbSNP_GRCh38.vcf.gz"
        idx = self.dbsnp_path / "dbSNP_GRCh38.rsidx"
        Resolver.rsidx_search(rsids, vcf, idx, self.coords_by_rsid["GRCh38"])

    @staticmethod
    def rsidx_search(rsids, vcf, idx, vardict):
        with sqlite3.connect(idx) as db:
            for line in rsidx.search.search(rsids, db, vcf):
                values = line.strip().split("\t")
                pos = int(values[1])
                rsid = values[2]
                assert ";" not in rsid
                vardict[rsid] = pos
