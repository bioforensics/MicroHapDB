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

from dataclasses import dataclass
from liftover import get_lifter
import pandas as pd
from pathlib import Path


class InvalidMarkerNameError(ValueError):
    pass


@dataclass
class VariantList:
    refr: str
    chrom: str
    pos: list[int]


class BaseMarker:
    def __init__(self, name, index, xref=None, source=None):
        self.name = BaseMarker.check_name(name)
        self.index = index
        self.xref = xref
        self.source = source
        self.positions = dict(GRCh37=list(), GRCh38=list())

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

    def check_chrom(self):
        name_chrom = self.name[2:4]
        if name_chrom.startswith("0"):
            name_chrom = name_chrom[1:]
        chrom = self.chrom.replace("chr", "")
        if chrom != name_chrom:
            raise ValueError(f"{self.name}: chromosome mismatch, {name_chrom} vs. {chrom}")

    def resolve(self, resolver):
        raise NotImplementedError()

    @property
    def all_rsids_present(self):
        if self.varref is None:
            return False
        if len(self.varref) < self.base.numvars:
            return False
        sum37 = sum([1 for rsid in self.varref if rsid in self.index.coords_by_rsid["GRCh37"]])
        sum38 = sum([1 for rsid in self.varref if rsid in self.index.coords_by_rsid["GRCh38"]])
        return sum37 == sum38 == self.base.numvars


class ExplicitMarker(BaseMarker):
    def __init__(self, name, index, varlist, xref=None, varref=None, source=None):
        super(ExplicitMarker, self).__init__(name, index, xref=xref, source=source)
        self.varref = varref
        if self.all_rsids_present:
            raise ValueError("use ImplicitMarker class when a full complement of RSIDs is available")
        self.check_chrom()
        self.resolve(varlist)

    @property
    def refr(self):
        return self.varlist.refr

    @property
    def chrom(self):
        return self.varlist.chrom

    @property
    def positions(self):
        return self.varlist.pos

    def resolve(self, varlist):
        alt_refr = "GRCh38" if varlist.refr == "GRCh37" else "GRCh37"
        varmap = self.index.position_mapping[varlist.refr]
        self.alt_varlist = VariantList(alt_refr, varlist.chrom, [varmap[p] for p in varlist.pos])
        self.varlist = varlist




class BaseMarkerDefinition:
    @staticmethod
    def from_csv(csvpath, source=None):
        markers = pd.read_csv(csvpath)
        for n, row in markers.iterrows():
            yield BaseMarkerDefinition.from_row(row, source)

    @classmethod
    def from_row(cls, row, source=None):
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
            source=source,
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
        source=None,
    ):
        self.name = BaseMarkerDefinition.check_name(name)
        self.numvars = BaseMarkerDefinition.check_num_vars(numvars)
        self.xref = xref
        self.refr = refr
        self.chrom = chrom
        self.positions = positions
        self.varref = varref
        self.source = source
        self.validate_definition()
        self.check_chrom()

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
        dbsnp_path = Path(dbsnp_path)
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

    def liftover(self, positions, lifter):
        return [lifter[self.chrom][pos][0][1] for pos in positions]

    @property
    def labpi(self):
        return self.name[4:].split("-")[0]


class CompleteMarkerDefinition:
    def __init__(self, base, resolver):
        self.base = base
        self.resolver = resolver
        self.positions = dict(GRCh37=list(), GRCh38=list())
        self._resolve()

    def _resolve(self):
        if self.all_rsids_present:
            self.positions["GRCh37"] = [self.resolver.coords_by_rsid["GRCh37"][rsid] for rsid in self.base.varref]
            self.positions["GRCh38"] = [self.resolver.coords_by_rsid["GRCh38"][rsid] for rsid in self.base.varref]
        elif self.base.refr == "GRCh37":
            self.positions["GRCh37"] = self.base.positions
            lifter = get_lifter("hg19", "hg38")
            self.positions["GRCh38"] = [lifter[self.base.chrom][pos][0][1] for pos in self.base.positions]
        else:
            assert self.base.positions is not None, self.base.name
            self.positions["GRCh38"] = self.base.positions
            lifter = get_lifter("hg38", "hg19")
            self.positions["GRCh37"] = [lifter[self.base.chrom][pos][0][1] for pos in self.base.positions]
        if self.base.varref is not None and len(self.base.varref) == self.base.numvars and not self.all_rsids_present:
            print(f"[{self.base.source.name}::{self.name}] all RSIDs provided, but some not found in dbSNP; falling back to explicitly provided positions")

    @property
    def all_rsids_present(self):
        if self.base.varref is None:
            return False
        if len(self.base.varref) < self.base.numvars:
            return False
        sum37 = sum([1 for rsid in self.base.varref if rsid in self.resolver.coords_by_rsid["GRCh37"]])
        sum38 = sum([1 for rsid in self.base.varref if rsid in self.resolver.coords_by_rsid["GRCh38"]])
        return sum37 == sum38 == self.base.numvars

    @property
    def region(self):
        if self.base.positions is None:
            return None
        start, end = min(self.base.positions), max(self.base.positions)
        return f"{self.base.refr}|{self.base.chrom}:{start}-{end}"

    @property
    def region37(self):
        return self._region(self.positions["GRCh37"])

    @property
    def region38(self):
        return self._region(self.positions["GRCh38"])

    def _region(self, poslist):
        start, end = min(poslist), max(poslist)
        return f"{self.base.chrom}:{start}-{end}"

    @property
    def is_match(self):
        if self.base.positions is None or not self.all_rsids_present:
            return None
        return sorted(self.base.positions) == sorted(self.positions[self.base.refr])

    def __str__(self):
        return f"{self.base.name}\t{self.is_match}\t{self.region}\t{self.region37}\t{self.region38}"

    @property
    def start(self):
        return min(self.positions["GRCh38"])

    @property
    def end(self):
        return max(self.positions["GRCh38"])

    @property
    def span(self):
        return (self.start, self.end)

    @property
    def posstr(self):
        return ";".join(map(str, self.positions["GRCh38"]))

    @property
    def name(self):
        return self.base.name

    @property
    def chrom_num(self):
        if self.base.chrom == "chrX":
            return 23
        else:
            return int(self.base.chrom[3:])
