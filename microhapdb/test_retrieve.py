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
    assert len(microhapdb.populations) == 96
    assert len(microhapdb.loci) == 198


def test_allele_frequencies():
    """Allele frequencies for 198 loci across 96 populations.

    >>> freqs = microhapdb.frequencies
    >>> alleles = freqs.query('Locus == "MHDBL000077"').Allele
    >>> list(alleles.unique())
    ['T,A,G', 'T,A,A', 'T,C,G', 'C,A,G']
    >>> pop_freqs = freqs.query('Locus == "MHDBL000077" and Allele == "T,A,A"')
    >>> len(pop_freqs)
    83
    >>> freqs.query('Locus == "MHDBL000077" and Allele == "T,A,A" and Population == "MHDBP000034"')
                 Locus   Population Allele  Frequency
    30638  MHDBL000077  MHDBP000034  T,A,A      0.107
    """
    af = microhapdb.frequencies
    assert af.shape == (82167, 4)
    result = af[af.Locus == 'MHDBL000135'].Allele.unique()
    assert len(result) == 6
    result = af[(af.Locus == 'MHDBL000135') & (af.Allele == 'G,G,G')]
    assert len(result) == 83
    result = af.query('Locus == "MHDBL000135" & Allele == "G,G,G" & Population == "MHDBP000092"').Frequency.values[0]
    assert pytest.approx(result, 0.202)


def test_loci():
    """Microhaplotype locus data

    >>> microhapdb.loci.query('ID == "MHDBL000077"')
                 ID Reference  Chrom     Start       End  Source
    76  MHDBL000077    GRCh38  chr15  63806358  63806495  ALFRED
    >>> microhapdb.id_xref('mh04CP-002')
                 ID Reference Chrom     Start       End  Source
    14  MHDBL000015    GRCh38  chr4  24304953  24304972  ALFRED
    """
    loc = microhapdb.loci
    vm = microhapdb.variantmap
    vr = microhapdb.variants
    assert loc.shape == (198, 6)
    result = loc[loc.Chrom == 'chr19']
    assert len(result) == 5
    varids = vm[vm.LocusID.isin(result.ID)]
    variants = vr[vr.ID.isin(varids.VariantID)]
    assert len(variants) == 17


def test_populations():
    """Population data

    >>> microhapdb.populations.query('ID == "MHDBP000048"')
                 ID     Name  Source
    47  MHDBP000048  Koreans  ALFRED
    >>> microhapdb.id_xref('SA004059S')
                 ID Name  Source
    73  MHDBP000074  Han  ALFRED
    >>> microhapdb.populations.query('Name.str.contains("Afr")')
                 ID               Name  Source
    40  MHDBP000041  African Americans  ALFRED
    67  MHDBP000068  African Americans  ALFRED
    83  MHDBP000084    Afro-Caribbeans  ALFRED
    """
    pop = microhapdb.populations
    assert pop.shape == (96, 3)
    assert pop[pop.ID == 'MHDBP000025'].Name.values == ['Karitiana']
    result = pop[pop.Name.str.contains('Jews')].ID.values
    assert list(result) == ['MHDBP000013', 'MHDBP000014', 'MHDBP000036', 'MHDBP000045']


def test_variants():
    """Microhaplotype variant data

    >>> microhapdb.id_xref('rs80047978')
                       ID Reference  Chrom  Position Alleles    Source
    13055  MHDBV000013056    GRCh38  chr15  63806494     A,G  dbSNP151
    >>> microhapdb.variants.query('Chrom == "chr15" and 52192400 < Position < 52192500')
                       ID Reference  Chrom  Position Alleles    Source
    12779  MHDBV000012780    GRCh38  chr15  52192466     C,T  dbSNP151
    12780  MHDBV000012781    GRCh38  chr15  52192467     G,T  dbSNP151
    12781  MHDBV000012782    GRCh38  chr15  52192471     C,T  dbSNP151
    12782  MHDBV000012783    GRCh38  chr15  52192490     C,T  dbSNP151
    12783  MHDBV000012784    GRCh38  chr15  52192495     C,T  dbSNP151
    """
    v = microhapdb.variants
    assert v.shape == (36304, 6)
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

    results = list(fetch_by_query('population', 'ID == "MHDBP000013"'))
    assert len(results) == 1
    assert results[0].Name.values == ['Jews, Ethiopian']


def test_fetch_by_id():
    assert list(fetch_by_id(None)) == []
    assert list(fetch_by_id('NotARealID')) == []
    results = list(fetch_by_id('rs1363241798'))
    assert len(results) == 1
    assert results[0].ID.values == ['MHDBV000022000']
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
    assert list(results[0].ID.values) == ['MHDBL000092']
    assert list(results[1].ID.values) == ['MHDBV000008871', 'MHDBV000008872',
                                          'MHDBV000008873', 'MHDBV000008874',
                                          'MHDBV000008875']
