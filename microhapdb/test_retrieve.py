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
    assert len(microhapdb.populations) == 96 + 3
    assert len(microhapdb.loci) == 198 + 15


def test_allele_frequencies():
    """Allele frequencies for 198 loci across 96 populations.

    >>> f = microhapdb.frequencies
    >>> alleles = f[f.Locus == 'MHDBL000077'].Allele
    >>> list(alleles.unique())
    ['A,G,A', 'A,G,C', 'A,A,C', 'C,G,C']
    >>> popfreq = f[(f.Locus == 'MHDBL000077') & (f.Allele == 'A,A,C')]
    >>> len(popfreq)
    26
    >>> f.query('Locus == "MHDBL000077" and Allele == "A,A,C" and Population == "MHDBP000025"')
                 Locus   Population Allele  Frequency
    50642  MHDBL000077  MHDBP000025  A,A,C      0.318
    """
    af = microhapdb.frequencies
    assert af.shape == (82533, 4)
    result = af[af.Locus == 'MHDBL000135'].Allele.unique()
    assert len(result) == 8
    result = af[(af.Locus == 'MHDBL000135') & (af.Allele == 'A,C,T')]
    assert len(result) == 96
    result = af.query('Locus == "MHDBL000135" & Allele == "A,C,T" & Population == "MHDBP000084"').Frequency.values[0]
    assert result == pytest.approx(0.025)


def test_loci():
    """Microhaplotype locus data

    >>> l = microhapdb.loci
    >>> l[l.ID == 'MHDBL000108']
                  ID Reference  Chrom     Start       End   AvgAe  Source
    107  MHDBL000108    GRCh38  chr18  78329886  78329968  3.2906  ALFRED
    >>> l[l.ID == 'MHDBL000065']
                 ID Reference  Chrom     Start       End   AvgAe Source
    64  MHDBL000065    GRCh38  chr14  32203269  32203324  3.7886   LOVD
    >>> microhapdb.id_xref('mh04CP-002')
                  ID Reference Chrom     Start       End   AvgAe  Source
    159  MHDBL000160    GRCh38  chr4  24304953  24304972  3.5182  ALFRED
    """
    loc = microhapdb.loci
    vm = microhapdb.variantmap
    vr = microhapdb.variants
    assert loc.shape == (213, 7)
    result = loc[loc.Chrom == 'chr19']
    assert len(result) == 5
    varids = vm[vm.LocusID.isin(result.ID)]
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
    assert pop.shape == (99, 3)
    assert pop[pop.ID == 'MHDBP000025'].Name.values == ['Finns']
    assert pop[pop.ID == 'MHDBP000049'].Name.values == ['Karitiana']
    result = pop[pop.Name.str.contains('Jews')].ID.values
    assert list(result) == ['MHDBP000044', 'MHDBP000045', 'MHDBP000046', 'MHDBP000047']


def test_variants():
    """Microhaplotype variant data

    >>> microhapdb.id_xref('SI338744D')
                      ID Reference Chrom   Position Alleles    Source
    1958  MHDBV000001959    GRCh38  chr1  161985865     A,G  dbSNP151
    >>> microhapdb.id_xref('rs80047978')
                       ID Reference  Chrom  Position Alleles    Source
    14355  MHDBV000014356    GRCh38  chr15  63806494     A,G  dbSNP151
    >>> microhapdb.variants.query('Chrom == "chr15" and 52192400 < Position < 52192500')
                       ID Reference  Chrom  Position Alleles    Source
    14079  MHDBV000014080    GRCh38  chr15  52192466     C,T  dbSNP151
    14080  MHDBV000014081    GRCh38  chr15  52192467     G,T  dbSNP151
    14081  MHDBV000014082    GRCh38  chr15  52192471     C,T  dbSNP151
    14082  MHDBV000014083    GRCh38  chr15  52192490     C,T  dbSNP151
    14083  MHDBV000014084    GRCh38  chr15  52192495     C,T  dbSNP151
    """
    v = microhapdb.variants
    assert v.shape == (40395, 6)
    assert len(v[v.Chrom == 'chr12']) == 1395


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
    assert results[0].ID.values == ['MHDBV000024096']
    assert results[0].Chrom.values == ['chr20']
    assert results[0].Position.values == [63539694]


def test_fetch_by_region():
    with pytest.raises(ValueError) as ve:
        list(fetch_by_region('chr5:1000000-2000000', 'population'))
    assert 'region query not supported for table "population"' in str(ve)

    with pytest.raises(ValueError) as ve:
        list(fetch_by_region('chr7:123-456-789', 'locus'))
    assert 'cannot parse region "chr7:123-456-789"' in str(ve)

    assert list(fetch_by_region('chrX', 'variant')) == []
    assert list(fetch_by_region('chrY', 'locus')) == []

    results = list(fetch_by_region('chr12:102866940-102866950', None))
    assert len(results) == 2
    print(results[0].ID.values)
    print(results[1].ID.values)
    assert list(results[0].ID.values) == ['MHDBL000053']
    assert list(results[1].ID.values) == ['MHDBV000009392', 'MHDBV000009393',
                                          'MHDBV000009394', 'MHDBV000009395',
                                          'MHDBV000009396']
