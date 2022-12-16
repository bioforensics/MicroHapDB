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

from liftover import get_lifter
import pandas as pd
from pathlib import Path
import rsidx
import sqlite3
from subprocess import run


class InvalidMarkerNameError(ValueError):
    pass


class Marker:
    @staticmethod
    def from_csv(csvpath, dbsnp_path=None):
        markers = pd.read_csv(csvpath)
        for n, row in markers.iterrows():
            yield Marker.from_row(row, dbsnp_path=dbsnp_path)
            print(".", flush=True, end="")

    @classmethod
    def from_row(cls, row, dbsnp_path=None):
        xref = None if pd.isna(row.Xref) else row.Xref
        refr = None if pd.isna(row.Refr) else row.Refr
        chrom = None if pd.isna(row.Chrom) else row.Chrom
        positions = None
        if not pd.isna(row.Positions):
            positions = str(row.Positions).split(";")
            positions = list(map(int, positions))
        rsids = None
        if not pd.isna(row.VarRef):
            rsids = row.VarRef.split(";")
        return cls(
            row.Name,
            row.NumVars,
            xref=xref,
            refr=refr,
            chrom=chrom,
            positions=positions,
            varref=rsids,
            dbsnp_path=dbsnp_path
        )

    def __init__(
        self,
        name,
        numvars,
        refr=None,
        chrom=None,
        positions=None,
        varref=None,
        xref=None,
        dbsnp_path=None,
    ):
        self.name = Marker.check_name(name)
        self.numvars = Marker.check_num_vars(numvars)
        self.xref = xref
        self.refr = refr
        self.chrom = chrom
        self.positions = positions
        self.varref = varref
        self.validate_definition()
        self.check_chrom()
        self.set_positions(dbsnp_path=Path(dbsnp_path))

    @staticmethod
    def check_name(name):
        strikes = list()
        if not name.startswith(("mh", "MH")):
            strikes.append("does not start with mh/MH")
        chromstr = name[2:4]
        try:
            chrom = int(chromstr)
            if chrom < 1 or chrom > 22:
                strikes.append(f"invalid chromosome '{chromstr}'")
        except:
            if chromstr != "0X":
                strikes.append(f"invalid chromosome '{chromstr}'")
        labpi, mhnumber = name[4:].split("-", 1)
        if len(labpi) < 2 or len(labpi) > 4:
            strikes.append(f"lab or PI designation '{labpi}' isn't between 2-4 characters")
        if "." in mhnumber:
            if mhnumber.count(".") > 1:
                strikes.append(f"invalid marker number '{mhnumber}'")
            version = mhnumber.split(".")[-1]
            if not version.startswith("v"):
                strikes.append(f"invalid version designator in marker number '{mhnumber}'")
            try:
                int(version[1:])
            except:
                strikes.append(f"invalid version designator in  marker number '{mhnumber}'")
        if len(strikes) > 0:
            message = f"{name}: " + "; ".join(strikes)
            raise InvalidMarkerNameError(message)
        return name.replace("MH", "mh")

    @staticmethod
    def check_num_vars(numvars):
        numvars = int(numvars)
        if numvars < 1:
            raise ValueError(f"invalid number of variants '{numvars}'")
        return numvars

    def validate_definition(self):
        if self.positions is None and self.varref is None:
            message = "must specify chromosome positions and/or RSIDs for each variant"
            raise ValueError(f"{self.name}: message")
        if self.positions is not None and (self.chrom is None or self.refr is None):
            message = f"{self.name}: if chromosome positions are explicitly specified, the reference assembly and chromosome must also be specified"
            raise ValueError(message)
        if self.positions is not None and len(self.positions) < self.numvars:
            message = f"{self.name}: expected {self.numvars} SNP/indel variants, but found only {len(self.positions)} positions specified"
            raise ValueError(message)
        if self.varref is not None and len(self.varref) < self.numvars:
            if self.positions is None or len(self.positions) < self.numvars:
                message = f"{self.name}: if a full complement of RSIDs cannot be specified, a full complement of positions must be specified; expected {self.numvars} SNP/indel variants, but found only {len(self.positions)} positions specified"
                raise ValueError(message)
        if self.refr not in (None, "GRCh37", "GRCh38"):
            raise ValueError(f"{self.name}: unsupported human genome assembly {self.refr}")

    def check_chrom(self):
        if self.chrom is None:
            return
        name_chrom = self.name[2:4]
        if name_chrom.startswith("0"):
            name_chrom = name_chrom[1:]
        chrom = self.chrom.replace("chr", "")
        if chrom != name_chrom:
            message = f"{self.name}: chromosome designation in marker name does not match user-supplied chromsome '{chrom}'"
            raise ValueError(message)

    def set_positions(self, dbsnp_path=None):
        if dbsnp_path is None:
            return
        if self.varref is not None and len(self.varref) == self.numvars:
            self.set_positions_from_rsids(dbsnp_path)
        else:
            assert len(self.positions) == self.numvars
            if self.refr == "GRCh37":
                self.positions37 = self.positions
                self.positions38 = self.liftover(self.positions, get_lifter("hg19", "hg38"))
            elif self.refr == "GRCh38":
                self.positions38 = self.positions
                self.positions37 = self.liftover(self.positions, get_lifter("hg38", "hg19"))
            else:
                raise ValueError(f"unsupported reference '{self.refr}'")

    def set_positions_from_rsids(self, dbsnp_path):
        self.positions37 = Marker.rsidx_search(self.varref, dbsnp_path / "dbSNP_GRCh37.vcf.gz", dbsnp_path / "dbSNP_GRCh37.rsidx")
        self.positions38 = Marker.rsidx_search(self.varref, dbsnp_path / "dbSNP_GRCh38.vcf.gz", dbsnp_path / "dbSNP_GRCh38.rsidx")

    @staticmethod
    def rsidx_search(rsids, vcf, idx):
        positions = dict()
        with sqlite3.connect(idx) as db:
            for line in rsidx.search.search(rsids, db, vcf):
                values = line.strip().split("\t")
                pos = int(values[1])
                rsid = values[2]
                assert ";" not in rsid
                positions[rsid] = pos
        return positions

    def liftover(self, positions, lifter):
        return [lifter[self.chrom][pos][0][1] for pos in positions]

    @property
    def labpi(self):
        return self.name[4:].split("-")[0]

    @property
    def region(self):
        if self.positions is None:
            return None
        start, end = min(self.positions), max(self.positions)
        return f"{self.refr}|{self.chrom}:{start}-{end}"

    @property
    def region37(self):
        return self._region("positions37")

    @property
    def region38(self):
        return self._region("positions38")

    def _region(self, poslist):
        if not hasattr(self, poslist):
            return None
        positions = getattr(self, poslist)
        if isinstance(positions, dict):
            positions = list(positions.values())
        start, end = min(positions), max(positions)
        return f"{self.chrom}:{start}-{end}"
