#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


import microhapdb
from microhapdb.retrieve import fetch_by_id, fetch_by_query, fetch_by_region
import pytest


def test_assumptions():
    assert len(microhapdb.populations) == 96 + 3 + 1
    assert len(microhapdb.markers) == 198 + 15 + 40


def test_allele_frequencies():
    """Allele frequencies for 198 markers across 96 populations.

    >>> f = microhapdb.frequencies
    >>> alleles = f[f.Marker == 'MHDBM000091'].Allele
    >>> list(alleles.unique())
    ['A,G,A', 'A,G,C', 'A,A,C', 'C,G,C']
    >>> popfreq = f[(f.Marker == 'MHDBM000091') & (f.Allele == 'A,A,C')]
    >>> len(popfreq)
    26
    >>> f.query('Marker == "MHDBM000091" and Allele == "A,A,C" and Population == "MHDBP000025"')
                Marker   Population Allele  Frequency
    50642  MHDBM000091  MHDBP000025  A,A,C      0.318
    """
    af = microhapdb.frequencies
    assert af.shape == (82670, 4)
    result = af[af.Marker == 'MHDBM000160'].Allele.unique()
    assert len(result) == 8
    result = af[(af.Marker == 'MHDBM000160') & (af.Allele == 'A,C,T')]
    assert len(result) == 96
    result = af.query('Marker == "MHDBM000160" & Allele == "A,C,T" & Population == "MHDBP000084"').Frequency.values[0]
    assert result == pytest.approx(0.025)


def test_markers():
    """Microhaplotype marker data

    >>> m = microhapdb.markers
    >>> m[m.ID == 'MHDBM000128']
                  ID Reference  Chrom     Start       End   AvgAe  Source
    127  MHDBM000128    GRCh38  chr18  78329885  78329968  3.2906  ALFRED
    >>> m[m.ID == 'MHDBM000078']
                 ID Reference  Chrom     Start       End   AvgAe                            Source
    77  MHDBM000078    GRCh38  chr14  32203273  32203324  3.7886  doi:10.1016/j.fsigen.2018.05.008
    >>> microhapdb.id_xref('mh04CP-002')
                  ID Reference Chrom     Start       End   AvgAe  Source
    192  MHDBM000193    GRCh38  chr4  24304952  24304972  3.5182  ALFRED
    """
    m = microhapdb.markers
    vm = microhapdb.variantmap
    vr = microhapdb.variants
    assert m.shape == (253, 7)
    result = m[m.Chrom == 'chr19']
    assert len(result) == 6
    varids = vm[vm.MarkerID.isin(result.ID)]
    variants = vr[vr.ID.isin(varids.VariantID)]
    assert len(variants) == 17


def test_populations():
    """Population data

    >>> p = microhapdb.populations
    >>> p[p.ID == 'MHDBP000048']
                 ID     Name  Source
    47  MHDBP000048  Kachari  ALFRED
    >>> p.query('ID == "MHDBP000053"')
                 ID     Name  Source
    52  MHDBP000053  Koreans  ALFRED
    >>> microhapdb.id_xref('SA004059S')
                 ID Name  Source
    32  MHDBP000033  Han  ALFRED
    >>> p.query('Name.str.contains("Afr")')
                ID               Name  Source
    1  MHDBP000002             Africa    LOVD
    2  MHDBP000003  African Americans  ALFRED
    3  MHDBP000004  African Americans  ALFRED
    4  MHDBP000005    Afro-Caribbeans  ALFRED
    """
    pop = microhapdb.populations
    assert pop.shape == (100, 3)
    assert pop[pop.ID == 'MHDBP000025'].Name.values == ['Finns']
    assert pop[pop.ID == 'MHDBP000049'].Name.values == ['Karitiana']
    result = pop[pop.Name.str.contains('Jews')].ID.values
    assert list(result) == ['MHDBP000044', 'MHDBP000045', 'MHDBP000046', 'MHDBP000047']


def test_variants():
    """Microhaplotype variant data

    >>> microhapdb.id_xref('SI338744D')
                      ID Reference Chrom   Position Alleles    Source
    1954  MHDBV000001955    GRCh38  chr1  161985865     A,G  dbSNP151
    >>> microhapdb.id_xref('rs80047978')
                       ID Reference  Chrom  Position Alleles    Source
    14094  MHDBV000014095    GRCh38  chr15  63806494     A,G  dbSNP151
    >>> microhapdb.variants.query('Chrom == "chr15" and 52192400 < Position < 52192500')
                       ID Reference  Chrom  Position Alleles    Source
    13818  MHDBV000013819    GRCh38  chr15  52192466     C,T  dbSNP151
    13819  MHDBV000013820    GRCh38  chr15  52192467     G,T  dbSNP151
    13820  MHDBV000013821    GRCh38  chr15  52192471     C,T  dbSNP151
    13821  MHDBV000013822    GRCh38  chr15  52192490     C,T  dbSNP151
    13822  MHDBV000013823    GRCh38  chr15  52192495     C,T  dbSNP151
    """
    v = microhapdb.variants
    assert v.shape == (39005, 6)
    assert len(v[v.Chrom == 'chr12']) == 1387


def test_fetch_by_query():
    with pytest.raises(ValueError) as ve:
        next(fetch_by_query(None, 'ID == "MHDBP000013"'))
    assert 'must specify table to invoke a query' in str(ve.value)

    with pytest.raises(ValueError) as ve:
        next(fetch_by_query('BoogerSnot', 'ID == "MHDBP000013"'))
    assert 'unsupported table "BoogerSnot"' in str(ve.value)

    results = list(fetch_by_query('population', 'ID == "BogusID"'))
    assert results == []

    results = list(fetch_by_query('population', 'ID == "MHDBP000045"'))
    assert len(results) == 1
    assert results[0].Name.values == ['Jews, Ethiopian']


def test_fetch_by_id():
    assert list(fetch_by_id(None)) == []
    assert list(fetch_by_id('NotARealID')) == []
    results = list(fetch_by_id('rs1363241798'))
    assert len(results) == 1
    print(results)
    assert results[0].ID.values == ['MHDBV000023609']
    assert results[0].Chrom.values == ['chr20']
    assert results[0].Position.values == [63539694]


def test_fetch_by_region():
    message = r'region query not supported for table "population"'
    with pytest.raises(ValueError, match=message) as ve:
        list(fetch_by_region('chr5:1000000-2000000', 'population'))

    with pytest.raises(ValueError, match='cannot parse region "chr7:123-456-789"') as ve:
        list(fetch_by_region('chr7:123-456-789', 'marker'))

    assert list(fetch_by_region('chrX', 'variant')) == []
    assert list(fetch_by_region('chrY', 'marker')) == []

    results = list(fetch_by_region('chr12:102866940-102866950', None))
    assert len(results) == 2
    print(results[0].ID.values)
    print(results[1].ID.values)
    assert list(results[0].ID.values) == ['MHDBM000062']
    assert list(results[1].ID.values) == [
        'MHDBV000009192', 'MHDBV000009193', 'MHDBV000009194', 'MHDBV000009195',
        'MHDBV000009196',
    ]


@pytest.mark.parametrize('markeraccession', [
    'mh13KK-218',
    'SI664607D',
    'MHDBM000072',
])
def test_allele_positions(markeraccession):
    pos = microhapdb.retrieve.allele_positions(markeraccession)
    assert pos == [53486691, 53486745, 53486756, 53486836]
    with pytest.raises(StopIteration):
        _ = microhapdb.retrieve.allele_positions('b0GusId')


def test_standardize_ids():
    from microhapdb.retrieve import standardize_ids as sid
    assert sid(['BoGUSid']) == list()
    assert sid(['MHDBM000114']) == ['MHDBM000114']
    assert sid(['SI664623B']) == ['MHDBM000114']
    assert sid(['MHDBM000114', 'SI664623B']) == ['MHDBM000114']
    assert sid(['SI664623B', 'NotARealId']) == ['MHDBM000114']
    assert sid(['SI664623B', 'mh04KK-011', 'MHDBM000135']) == ['MHDBM000114', 'MHDBM000135', 'MHDBM000194']
    assert sid(['rs547950691', 'mh02KK-131', 'SA002765U']) == ['MHDBM000140', 'MHDBP000054', 'MHDBV000027667']
