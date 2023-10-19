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
import pytest


def test_assumptions():
    assert microhapdb.repeats.shape == (1380, 1)


@pytest.mark.parametrize(
    "marker,repeat",
    [
        ("mh22USC-22qB.v2", "L1ME3A|LINE|L1"),
        ("mh08SCUZJ-0484105", "AluSx3|SINE|Alu;AluYa5|SINE|Alu"),
        ("mh13WL-037.v1", "L2|LINE"),
        ("mh06SCUZJ-0001113", "MER44C|DNA|TcMar-Tigger"),
    ],
)
def test_repeats_basic(marker, repeat):
    assert microhapdb.repeats.loc[marker].Repeat == repeat
