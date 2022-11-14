# -------------------------------------------------------------------------------------------------
# Copyright (c) 2018, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------


from microhapdb import retrieve_by_id
import pytest


def test_retrieve_by_id():
    with pytest.raises(ValueError, match=r'identifier "NotARealID" not found in MicroHapDB'):
        retrieve_by_id("NotARealID")
    results = retrieve_by_id("rs62069897")
    assert len(results) == 1
    assert results.Name.values == ["mh17CP-006"]
    assert results.Chrom.values == ["chr17"]
    assert results.Offsets.values == ["15385768,15385781,15385817"]


def test_all_ids():
    with pytest.raises(ValueError, match=r'identifier "BoGUSid" not found in MicroHapDB'):
        retrieve_by_id("BoGUSid")
    assert retrieve_by_id("mh03KK-150").Name.tolist() == ["mh03KK-150"]
    assert retrieve_by_id("SI664623B").Name.tolist() == ["mh17KK-077"]
    assert retrieve_by_id("MHDBM-c2b1818c").Name.tolist() == ["mh15CP-003"]
    assert retrieve_by_id("MHDBM-f4474b7c").Name.tolist() == ["mh17KK-054", "mh17AT-35"]
    assert retrieve_by_id("rs58111155").Name.tolist() == ["mh01KKCS-001", "mh01KK-001"]
    assert retrieve_by_id("PEL").ID.tolist() == ["PEL"]
    assert retrieve_by_id("Han").ID.tolist() == ["MHDBP-48c2cfb2aa", "SA000001B", "SA000009J"]
