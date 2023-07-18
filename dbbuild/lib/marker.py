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

from .variant import VariantList
import pandas as pd


class InvalidMarkerNameError(ValueError):
    pass


class Marker:
    field_names = (
        "Name",
        "NumVars",
        "Extent",
        "Chrom",
        "Start",
        "End",
        "Positions",
        "Positions37",
        "RSIDs",
        "Source",
    )

    @staticmethod
    def from_csv(csvpath, index, source=None):
        table = pd.read_csv(csvpath)
        for n, row in table.iterrows():
            rsids = []
            if not pd.isna(row.VarRef):
                rsids = row.VarRef.split(";")
            if len(rsids) == row.NumVars:
                yield MarkerFromIDs(
                    row.Name, row.Chrom, rsids, index, xrefs=row.Xref, source=source
                )
            else:
                positions = row.Positions.split(";")
                positions = list(map(int, positions))
                position_list = VariantList(row.Refr, row.Chrom, positions)
                yield MarkerFromPositions(
                    row.Name, position_list, rsids, index, xrefs=row.Xref, source=source
                )

    def __init__(self, name, rsids, index, xrefs=None, source=None):
        self.name = Marker.check_name(name)
        self.rsids = rsids
        self.index = index
        self.xrefs = xrefs
        self.sources = list()
        if source is not None:
            self.sources.append(source)
        self.positions = dict(GRCh37=list(), GRCh38=list())

    def __len__(self):
        return self.end - self.start + 1

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
        if len(labpi) < 2 or len(labpi) > 6:
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
                strikes.append(f"invalid version designator in marker number '{mhnumber}'")
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
    def numvars(self):
        raise NotImplementedError()

    @property
    def all_rsids_present(self):
        if self.rsids is None:
            return False
        if len(self.rsids) < self.numvars:
            return False
        num37 = len(list(self.index.resolve(self.rsids, "GRCh37")))
        num38 = len(list(self.index.resolve(self.rsids, "GRCh38")))
        return num37 == num38 == self.base.numvars

    @property
    def fields(self):
        pos37 = self.posstr("GRCh37")
        pos38 = self.posstr("GRCh38")
        rsids = ";".join(self.rsids)
        return (
            self.name,
            self.numvars,
            len(self),
            self.chrom,
            self.start,
            self.end,
            pos38,
            pos37,
            rsids,
            self.sourcename,
        )

    @property
    def start(self):
        return min(self.positions["GRCh38"])

    @property
    def end(self):
        return max(self.positions["GRCh38"])

    @property
    def span(self):
        return self.start, self.end

    @property
    def sourcename(self):
        if len(self.sources) == 0:
            return None
        names = [s.name for s in sorted(self.sources, key=lambda s: s.sortkey)]
        return ";".join(names)

    @property
    def source(self):
        return self.sources[0]

    def posstr(self, refr="GRCh38"):
        return ";".join(map(str, self.positions[refr]))

    @property
    def chrom_num(self):
        if self.chrom == "chrX":
            return 23
        return int(self.chrom[3:])

    @property
    def locus(self):
        return self.name.split(".")[0]

    def overlaps(self, other):
        same_chrom = self.chrom_num == other.chrom_num
        return same_chrom and self.start <= other.end and self.end >= other.start

    @property
    def sortkey(self):
        return self.chrom_num, self.span, self.name


class MarkerFromPositions(Marker):
    def __init__(self, name, positions, rsids, index, xrefs=None, source=None):
        super(MarkerFromPositions, self).__init__(name, rsids, index, xrefs=xrefs, source=source)
        self.chrom = positions.chrom
        self.check_chrom()
        self.resolve(positions)
        if self.all_rsids_present:
            raise ValueError(
                "use MarkerFromIDs class when a full complement of RSIDs is available"
            )

    def resolve(self, varlist):
        self.positions[varlist.refr] = varlist.pos
        alt_refr = "GRCh38" if varlist.refr == "GRCh37" else "GRCh37"
        self.positions[alt_refr] = list(
            self.index.map(varlist.refr, varlist.chrom, varlist.pos, strict=True)
        )

    @property
    def numvars(self):
        return len(self.positions["GRCh38"])


class MarkerFromIDs(Marker):
    def __init__(self, name, chrom, rsids, index, xrefs=None, source=None):
        super(MarkerFromIDs, self).__init__(name, rsids, index, xrefs=xrefs, source=source)
        self.chrom = chrom
        self.check_chrom()
        self.resolve()

    def resolve(self):
        for refr in ("GRCh37", "GRCh38"):
            for position in self.index.resolve(self.rsids, refr, strict=True):
                self.positions[refr].append(position)
            self.positions[refr].sort()

    @property
    def numvars(self):
        return len(self.rsids)
