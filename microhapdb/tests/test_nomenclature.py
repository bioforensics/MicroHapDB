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


import microhapdb
from microhapdb.nomenclature import Identifier, ValidationErrors, ValidationWarnings, validate
import pytest


def test_prefix():
    id1 = Identifier("mh01KK-001")
    assert id1.valid
    assert id1.chrom == "chr1"
    assert id1.warnings is None
    id2 = Identifier("MH21KK-320.v1")
    assert id2.valid
    assert ValidationWarnings.MHCASE in id2._validation_warnings
    id3 = Identifier("05WL-064")
    assert not id3.valid
    assert ValidationErrors.PREFIX in id3._validation_errors


def test_chrom():
    id1 = Identifier("mh10CP-006")
    assert id1.valid
    assert id1.chrom == "chr10"
    assert id1.warnings is None
    id2 = Identifier("mh25DSS-004")
    assert not id2.valid
    assert ValidationErrors.CHROM in id2._validation_errors
    id3 = Identifier("mhZZKK-123.v1")
    assert not id3.valid
    assert ValidationErrors.CHROM in id3._validation_errors
    id4 = Identifier("mh0XUSC-XpB")
    assert id4.valid
    assert id4.chrom == "chrX"
    assert id4.chrom_label == "0X"
    assert id4.chrom_num == 23
    id5 = Identifier("mh0YDSS-bogus")
    assert id5.valid
    assert id5.chrom == "chrY"
    assert id5.chrom_label == "0Y"
    assert id5.chrom_num == 24


def test_hyphen():
    id1 = Identifier("mh07KK009.v2")
    assert not id1.valid
    assert ValidationErrors.HYPHEN in id1._validation_errors
    id2 = Identifier("mh19-ZHA-006")
    assert not id2.valid
    assert ValidationErrors.HYPHEN in id2._validation_errors
    id3 = Identifier("mh04WL-079")
    assert id3.valid


def test_lab():
    id1 = Identifier("mh16DSSAGS-123")
    assert not id1.valid
    assert ValidationErrors.LABLENGTH in id1._validation_errors
    id2 = Identifier("mh13dss-xyz")
    assert id2.valid
    assert ValidationWarnings.LABCASE in id2._validation_warnings
    id3 = Identifier("mh21WL-21:20110:5")
    assert not id3.valid
    assert ValidationErrors.ALPHANUM in id3._validation_errors
    id4 = Identifier("mh04WL-098")
    assert id4.valid
    assert id4.warnings is None


def test_suffix():
    id1 = Identifier("mh06USC-6qD")
    assert id1.valid
    assert id1.suffix is None
    id2 = Identifier("mh16WL-003.v1")
    assert id2.valid
    assert id2.suffix == "v1"
    assert id2.locus == "mh16WL-003"
    id3 = Identifier("mh10NH-14.v1.4")
    assert not id3.valid
    assert ValidationErrors.PERIODS in id3._validation_errors
    id4 = Identifier("mh01KK-212.2")
    assert not id4.valid
    assert ValidationErrors.SUFFIXV in id4._validation_errors
    id5 = Identifier("mh08ZHA-011.xyz")
    assert not id5.valid
    assert ValidationErrors.SUFFIXNUM in id5._validation_errors


def test_warnings_and_errors():
    id1 = Identifier("MH19dss-abc")
    assert id1.valid
    assert id1.warnings == "'mh' prefix is not lower case; lab symbol is not upper case"
    id2 = Identifier("mh00DSS-007-5.v4")
    assert not id2.valid
    assert id2.errors == "invalid chromosome; too many or too few hypens"


def test_chrom_num():
    ids = [Identifier(idstr) for idstr in ("mh06KK-101", "mh02USC-2qB.v2", "mh17WL-036")]
    ids.sort(key=lambda x: x.chrom_num)
    assert [x.chrom_num for x in ids] == [2, 6, 17]


def test_invalid_attributes():
    mhid = Identifier("GATTACA")
    assert not mhid.valid
    assert mhid.chrom is None
    assert mhid.chrom_num is None
    assert mhid.chrom_label is None
    assert mhid.lab is None
    assert mhid.designation is None
    assert mhid.suffix is None
    assert mhid.locus is None
    assert str(mhid) == "invalid"
    assert mhid.detail is None


def test_validate_good(capsys):
    assert validate("mh13ZBF-001") is True
    terminal = capsys.readouterr()
    observed = terminal.out
    expected = """
mh13ZBF-001
├── mh  [standard prefix]
├── 13  [chromosome]
├── ZBF [lab]
├── -   [standard separator]
└── 001 [designation]
"""
    assert observed.strip() == expected.strip()
    assert terminal.err == ""


def test_validate_good_warning(capsys):
    assert validate("MH21KK-324.v2") is True
    terminal = capsys.readouterr()
    observed = terminal.out
    expected = """
mh21KK-324.v2
├── mh  [standard prefix]
├── 21  [chromosome]
├── KK  [lab]
├── -   [standard separator]
├── 324 [designation]
└── v2  [marker suffix]
    """
    assert observed.strip() == expected.strip()
    assert "[nomenclature] validation warnings: 'mh' prefix is not lower case" in terminal.err


def test_validate_bad(capsys):
    assert validate("MH35DSS-999.v4.5") is False
    terminal = capsys.readouterr()
    message = "[nomenclature] validation errors: invalid chromosome; contains too many period/full stop characters"
    assert message in terminal.err
    assert "[nomenclature] validation warnings: 'mh' prefix is not lower case" in terminal.err


@pytest.mark.parametrize("name", microhapdb.markers.Name)
def test_valid(name):
    ident = Identifier(name)
    print(ident.errors)
    assert ident.valid
