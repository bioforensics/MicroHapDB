# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


import microhapdb
from microhapdb.retrieve import by_id, by_region
from microhapdb.retrieve import standardize_marker_ids, standardize_population_ids
import pytest


def test_retrieve_by_id():
    with pytest.raises(ValueError, match=r'identifier "NotARealID" not found in MicroHapDB'):
        by_id('NotARealID')
    results = by_id('rs62069897')
    assert len(results) == 1
    print(results)
    assert results.Name.values == ['mh17CP-006']
    assert results.Chrom.values == ['chr17']
    assert results.Offsets.values == ['15385768,15385781,15385817']


def test_all_ids():
    with pytest.raises(ValueError, match=r'identifier "BoGUSid" not found in MicroHapDB'):
        by_id('BoGUSid')
    assert by_id('mh03KK-150').Name.values == ['mh03KK-150']
    assert by_id('SI664623B').Name.values == ['mh17KK-077']
    assert by_id('MHDBM-c2b1818c').Name.values == ['mh15CP-003']
    assert list(by_id('MHDBM-f4474b7c').Name.values) == ['mh17KK-054', 'mh17AT-35']
    assert by_id('rs58111155').Name.values == ['mh01KK-001']
    assert by_id('Peruvian').ID.values == ['SA004245P']
    assert list(by_id('Han').ID.values) == ['SA004058R', 'SA004059S', 'MHDBP-48c2cfb2aa', 'SA000009J', 'SA000001B']


def test_retrieve_by_region():
    with pytest.raises(ValueError, match='cannot parse region "chr7:123-456-789"'):
        list(by_region('chr7:123-456-789'))

    assert len(by_region('chrX')) == 0
    assert len(by_region('chrY')) == 0

    results = by_region('chr12:100000000-200000000')
    assert len(results) == 6
    print(results.Name.values)
    assert sorted(results.Name.values) == sorted([
        'mh12KK-093', 'mh12KK-045', 'mh12KK-042', 'mh12KK-046', 'mh12AT-25', 'mh12CP-003'
    ])


def test_marker_ids():
    assert len(standardize_marker_ids(['BoGUSid'])) == 0
    assert list(standardize_marker_ids(['mh05KK-058'])) == ['mh05KK-058']
    assert list(standardize_marker_ids(['SI664549I'])) == ['mh01KK-117']
    assert list(standardize_marker_ids(['MHDBM-3d69621c'])) == ['mh11KK-040', 'mh11AT-23']
