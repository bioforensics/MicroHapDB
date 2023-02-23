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
    assert results.Positions.values == ["15385769;15385782;15385818"]


def test_all_ids():
    with pytest.raises(ValueError, match=r'identifier "BoGUSid" not found in MicroHapDB'):
        retrieve_by_id("BoGUSid")
    assert retrieve_by_id("mh03KK-150.v1").Name.tolist() == ["mh03KK-150.v1"]
    assert retrieve_by_id("rs58111155").Name.tolist() == [
        "mh01KK-001.v3",
        "mh01KK-001.v5",
        "mh01KK-001.v1",
        "mh01KK-001.v4",
    ]
    assert retrieve_by_id("EUR").ID.tolist() == ["EUR"]
    assert sorted(retrieve_by_id("Han").ID.tolist()) == [
        "MHDBP-48c2cfb2aa",
        "SA000001B",
        "SA000009J",
    ]
