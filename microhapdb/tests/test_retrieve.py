# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


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
