# -------------------------------------------------------------------------------------------------
# Copyright (c) 2023, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from enum import Enum
from io import StringIO
import sys


class ValidationErrors(Enum):
    PREFIX = "required 'mh' prefix missing"
    CHROM = "invalid chromosome"
    HYPHEN = "too many or too few hypens"
    LABLENGTH = "lab symbol is not between 2-4 characters in length"
    ALPHANUM = "designation is not alphanumeric"
    PERIODS = "contains too many period/full stop characters"
    SUFFIXV = "suffix does not begin with 'v'"
    SUFFIXNUM = "suffix is not numeric"


class ValidationWarnings(Enum):
    MHCASE = "'mh' prefix is not lower case"
    LABCASE = "lab symbol is not upper case"


class Identifier:
    """Class for parsing and validating MH identifiers

    >>> mhid = Identifier("NotARealID")
    >>> mhid.valid
    False
    >>> mhid.errors
    "required 'mh' prefix missing; invalid chromosome; too many or too few hypens"
    >>> mhid.warnings
    >>> mhid = Identifier("mh14WL-003")
    >>> mhid.valid
    True
    >>> mhid.chrom
    'chr14'
    >>> mhid.lab
    'WL'
    >>> mhid.designation
    '003'
    >>> mhid.suffix
    >>> str(mhid)
    'mh14WL-003'
    >>> print(mhid.detail)
    mh14WL-003
    ├── mh  [standard prefix]
    ├── 14  [chromosome]
    ├── WL  [lab]
    ├── -   [standard separator]
    └── 003 [designation]
    >>> mhid = Identifier("mh12KK-202.v2")
    >>> mhid.lab
    'KK'
    >>> mhid.suffix
    'v2'
    >>> print(mhid.detail)
    mh12KK-202.v2
    ├── mh  [standard prefix]
    ├── 12  [chromosome]
    ├── KK  [lab]
    ├── -   [standard separator]
    ├── 202 [designation]
    └── v2  [marker suffix]
    """

    def __init__(self, mhid):
        self._raw = mhid
        self._validation_errors = list()
        self._validation_warnings = list()
        self._valid = self.validate()

    def validate(self):
        self.validate_prefix()
        if self._raw.count("-") != 1:
            self._validation_errors.append(ValidationErrors.HYPHEN)
            return False
        self.validate_lab()
        self.validate_suffix()
        return len(self._validation_errors) == 0

    def validate_prefix(self):
        prefix = self._raw[0:2]
        if prefix in ("mh", "MH"):
            if prefix == "MH":
                self._validation_warnings.append(ValidationWarnings.MHCASE)
        else:
            self._validation_errors.append(ValidationErrors.PREFIX)
        chrom = self._raw[2:4]
        try:
            chrom_num = int(chrom)
            if chrom_num < 1 or chrom_num > 22:
                self._validation_errors.append(ValidationErrors.CHROM)
        except:
            if chrom not in ("0X", "0Y"):
                self._validation_errors.append(ValidationErrors.CHROM)

    def validate_lab(self):
        lab = self._raw[4:].split("-")[0]
        if not lab.isupper():
            self._validation_warnings.append(ValidationWarnings.LABCASE)
        if len(lab) < 2 or len(lab) > 6:
            self._validation_errors.append(ValidationErrors.LABLENGTH)
        designation = self._raw.split("-")[1].split(".")[0]
        if not designation.isalnum():
            self._validation_errors.append(ValidationErrors.ALPHANUM)

    def validate_suffix(self):
        periodcount = self._raw.count(".")
        if periodcount > 1:
            self._validation_errors.append(ValidationErrors.PERIODS)
        elif periodcount == 1:
            suffix = self._raw.split(".")[1]
            if suffix[0] != "v":
                self._validation_errors.append(ValidationErrors.SUFFIXV)
            if len(suffix) > 1 and not suffix[1:].isdigit():
                self._validation_errors.append(ValidationErrors.SUFFIXNUM)

    @property
    def valid(self):
        return self._valid

    @property
    def errors(self):
        if self.valid:
            return None
        return "; ".join([e.value for e in self._validation_errors])

    @property
    def warnings(self):
        if len(self._validation_warnings) == 0:
            return None
        return "; ".join([w.value for w in self._validation_warnings])

    @property
    def chrom(self):
        if not self.valid:
            return None
        chrom = self._raw[2:4]
        if chrom[0] == "0":
            chrom = self._raw[3]
        return f"chr{chrom}"

    @property
    def chrom_num(self):
        if not self.valid:
            return None
        chrom = self._raw[2:4]
        if chrom == "0X":
            return 23
        elif chrom == "0Y":
            return 24
        else:
            return int(chrom)

    @property
    def chrom_label(self):
        if not self.valid:
            return None
        return self._raw[2:4]

    @property
    def lab(self):
        if not self.valid:
            return None
        return self._raw[4:].split("-")[0].upper()

    @property
    def designation(self):
        if not self.valid:
            return None
        return self._raw.split("-")[1].split(".")[0]

    @property
    def suffix(self):
        if not self.valid:
            return None
        if "." not in self._raw:
            return None
        return self._raw.split(".")[1]

    @property
    def locus(self):
        if not self.valid:
            return None
        return self._raw.split(".")[0]

    def __str__(self):
        if not self.valid:
            return "invalid"
        ident = f"mh{self.chrom_label}{self.lab}-{self.designation}"
        if self.suffix is not None:
            ident += f".{self.suffix}"
        return ident

    @property
    def detail(self):
        if not self.valid:
            return None
        lengths = [len(value) for value in (self.chrom_label, self.lab, self.designation)]
        if self.suffix is not None:
            lengths.append(len(self.suffix))
        length = max(lengths)
        strfmt = f"{{:{length}s}}"
        output = StringIO()
        print(str(self), file=output)
        print("├──", strfmt.format("mh"), "[standard prefix]", file=output)
        print("├──", strfmt.format(self.chrom_label), "[chromosome]", file=output)
        print("├──", strfmt.format(self.lab), "[lab]", file=output)
        print("├──", strfmt.format("-"), "[standard separator]", file=output)
        if self.suffix is None:
            print("└──", strfmt.format(self.designation), "[designation]", file=output)
        else:
            print("├──", strfmt.format(self.designation), "[designation]", file=output)
            print("└──", strfmt.format(self.suffix), "[marker suffix]", file=output)
        return output.getvalue()


def validate(ident):
    mhid = Identifier(str(ident))
    if mhid.valid:
        print(mhid.detail)
    else:
        print("[nomenclature] validation errors:", mhid.errors, file=sys.stderr)
    if mhid.warnings is not None:
        print("[nomenclature] validation warnings:", mhid.warnings, file=sys.stderr)
    return mhid.valid
