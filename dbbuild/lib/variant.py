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

from collections import defaultdict
from dataclasses import dataclass
import pandas as pd
from pathlib import Path
import rsidx
import sqlite3
import subprocess
from tempfile import TemporaryDirectory
from warnings import warn


@dataclass
class VariantList:
    refr: str
    chrom: str
    pos: list[int]

    def __len__(self):
        return len(self.pos)

    def __str__(self):
        return f"{self.refr}::{self.chrom}:{min(self.pos)}-{max(self.pos)}"


class VariantIndex:
    def __init__(self, table, dbsnp_path, chain_path):
        self.table = table
        self.dbsnp_path = Path(dbsnp_path)
        self.chain_path = Path(chain_path)
        self.coords_by_rsid = dict(GRCh37=dict(), GRCh38=dict())
        self.position_mapping = dict(GRCh37=defaultdict(dict), GRCh38=defaultdict(dict))
        self.resolve_all_rsids()
        self.map_all_positions()

    @staticmethod
    def table_from_filenames(filenames):
        return pd.concat([pd.read_csv(fn) for fn in filenames]).reset_index()

    def resolve_all_rsids(self):
        self.load_merged_rsids()
        rsids = set(self.all_rsids())
        vcf = self.dbsnp_path / "dbSNP_GRCh37.vcf.gz"
        idx = self.dbsnp_path / "dbSNP_GRCh37.rsidx"
        VariantIndex.rsidx_search(rsids, vcf, idx, self.coords_by_rsid["GRCh37"])
        vcf = self.dbsnp_path / "dbSNP_GRCh38.vcf.gz"
        idx = self.dbsnp_path / "dbSNP_GRCh38.rsidx"
        VariantIndex.rsidx_search(rsids, vcf, idx, self.coords_by_rsid["GRCh38"])

    def all_rsids(self):
        for n, row in self.table.iterrows():
            if pd.isna(row.VarRef):
                continue
            for rsid in row.VarRef.split(";"):
                yield rsid
                if rsid in self.merged_rsids:
                    yield self.merged_rsids[rsid]

    @staticmethod
    def rsidx_search(rsids, vcf, idx, vardict):
        if not idx.is_file():
            raise FileNotFoundError(idx)
        with sqlite3.connect(idx) as db:
            for line in rsidx.search.search(rsids, db, vcf):
                values = line.strip().split("\t")
                pos = int(values[1])
                rsid = values[2]
                assert ";" not in rsid
                vardict[rsid] = pos

    def load_merged_rsids(self, updateint=1e6):
        merged_file = self.dbsnp_path / "refsnp-merged.csv.gz"
        if not merged_file.is_file():
            merged_file = self.dbsnp_path / "refsnp-merged.csv"
        if merged_file:
            table = pd.read_csv(merged_file)
            self.merged_rsids = dict(zip(table.Source, table.Target))
        else:
            merged_file = self.dbsnp_path / "refsnp-merged.json"
            if not merged_file.is_file():
                raise FileNotFoundError(merged_file)
            self.merged_rsids = dict()
            threshold = updateint
            for n, line in enumerate(instream):
                try:
                    data = json.loads(line)
                except:
                    warn(f"Could not parse line {n+1}, skipping: {line}")
                source = data["refsnp_id"]
                targets = data["merged_snapshot_data"]["merged_into"]
                for target in targets:
                    self.merged_rsids[f"rs{source}"] = f"rs{target}"
                if n >= threshold:
                    threshold += updateint
                    if threshold == updateint * 10:
                        updateint = threshold
                    print(f"processed {n} rows")
            table = pd.DataFrame(self.merged_rsids.items(), columns=["Source", "Target"])
            table.to_csv(self.dbsnp_path / "refsnp-merged.csv", index=False)

    def map_all_positions(self):
        self.map_all_positions_for_refr("GRCh37", self.chain_path / "hg19ToHg38.over.chain.gz")
        self.map_all_positions_for_refr("GRCh38", self.chain_path / "hg38ToHg19.over.chain.gz")

    def map_all_positions_for_refr(self, refr, chain):
        with TemporaryDirectory() as tmpdir:
            source_file = Path(tmpdir) / f"{refr}-source.bed"
            dest_file = Path(tmpdir) / f"{refr}-dist.bed"
            unmapped_file = Path(tmpdir) / f"{refr}-unmapped.bed"
            source_positions = pd.DataFrame(
                self.all_positions(refr), columns=["Chrom", "Start", "End"]
            )
            source_positions.to_csv(source_file, sep="\t", header=False, index=False)
            args = map(str, ["liftOver", source_file, chain, dest_file, unmapped_file])
            subprocess.run(args)
            assert unmapped_file.stat().st_size == 0
            dest_positions = pd.read_csv(dest_file, sep="\t", names=["Chrom", "Start", "End"])
            for source, dest in zip(source_positions.iterrows(), dest_positions.iterrows()):
                self.position_mapping[refr][source[1].Chrom][source[1].End] = dest[1].End

    def all_positions(self, refr):
        for n, row in self.table.iterrows():
            if pd.isna(row.Refr) or row.Refr != refr:
                continue
            for position in row.Positions.split(";"):
                yield row.Chrom, int(position) - 1, int(position)

    def resolve(self, rsids, refr, strict=False):
        for rsid in rsids:
            while rsid not in self.coords_by_rsid[refr] and rsid in self.merged_rsids:
                print(f"swapping out {rsid} for {self.merged_rsids[rsid]}")
                rsid = self.merged_rsids[rsid]
            if rsid in self.coords_by_rsid[refr]:
                yield self.coords_by_rsid[refr][rsid]
            elif strict is True:
                raise ValueError(f"no {refr} coordinate for {rsid}")

    def map(self, refr, chrom, positions, strict=False):
        for position in positions:
            if position in self.position_mapping[refr][chrom]:
                yield self.position_mapping[refr][chrom][position]
            elif strict is True:
                raise ValueError(
                    f"no mapping for {chrom}:{position} from {refr} to alternate assembly"
                )
