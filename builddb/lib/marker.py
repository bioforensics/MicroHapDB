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


class InvalidMarkerNameError(ValueError):
    pass


class Marker:
    @staticmethod
    def from_csv(csvpath):
        markers = pd.read_csv(csvpath)
        for n, row in markers.iterrows():
            yield Marker.from_row(row)

    @classmethod
    def from_row(cls, row):
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
        )

    def __init__(
        self, name, numvars, refr=None, chrom=None, positions=None, varref=None, xref=None
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
        if positions is None and varref is None:
            message = f"{name}: must provide chromosome positions and/or RSIDs for each variant"
            raise ValueError(message)

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
        if self.refr not in (None, "GRCh37", "GRCh38", "CHM13"):
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

    @property
    def labpi(self):
        return self.name[4:].split("-")[0]
