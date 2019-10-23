# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


import microhapdb
from microhapdb.panel import panel_alpha, panel_alpha_legacy, panel_beta, panel_beta_legacy
import pytest


def test_alpha():
    assert len(panel_alpha_legacy()) == 22
    assert sorted(panel_alpha()) == [
        'mh01KK-117', 'mh02KK-134', 'mh03CP-005', 'mh04CP-002', 'mh05KK-170', 'mh06KK-008',
        'mh07CP-004', 'mh08KK-039', 'mh09KK-157', 'mh10KK-163', 'mh11KK-180', 'mh12CP-008',
        'mh13KK-218', 'mh14CP-003', 'mh15CP-001', 'mh16KK-049', 'mh17CP-001', 'mh18CP-005',
        'mh19KK-299', 'mh20KK-307', 'mh21KK-320', 'mh22KK-061'
    ]


def test_beta():
    assert len(panel_beta_legacy()) == 50
    assert sorted(panel_beta()) == [
        'mh01CP-016', 'mh01KK-117', 'mh01KK-205', 'mh02KK-136', 'mh02KK-138', 'mh03CP-005',
        'mh03KK-007', 'mh03KK-150', 'mh04KK-013', 'mh04KK-017', 'mh05KK-020', 'mh06CP-007',
        'mh06KK-008', 'mh07CP-004', 'mh07KK-030', 'mh07KK-031', 'mh08KK-032', 'mh08KK-039',
        'mh09KK-033', 'mh09KK-153', 'mh09KK-157', 'mh10CP-003', 'mh10KK-101', 'mh10KK-170',
        'mh11KK-037', 'mh11KK-180', 'mh11KK-191', 'mh12CP-008', 'mh12KK-046', 'mh12KK-202',
        'mh13KK-213', 'mh13KK-223', 'mh14CP-003', 'mh14KK-048', 'mh15KK-069', 'mh15KK-104',
        'mh16KK-049', 'mh16KK-255', 'mh17CP-001', 'mh17KK-052', 'mh17KK-055', 'mh18CP-005',
        'mh18KK-293', 'mh19KK-057', 'mh19KK-299', 'mh20KK-058', 'mh20KK-307', 'mh21KK-315',
        'mh22KK-060', 'mh22KK-061'
    ]
