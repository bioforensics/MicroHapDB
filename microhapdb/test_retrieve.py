#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


import microhapdb
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
    >>> microhapdb.fetch_by_id('mh04CP-002')
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
    >>> microhapdb.fetch_by_id('SA004059S')
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

    >>> microhapdb.fetch_by_id('rs80047978')
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


def test_query_mode(capsys):
    with pytest.raises(ValueError) as ve:
        microhapdb.retrieve.query_mode(None, 'ID == "MHDBP000013"')
    assert 'must specify table to invoke a query' in str(ve.value)

    with pytest.raises(ValueError) as ve:
        microhapdb.retrieve.query_mode('BoogerSnot', 'ID == "MHDBP000013"')
    assert 'unsupported table "BoogerSnot"' in str(ve.value)

    microhapdb.retrieve.query_mode('population', 'ID == "BogusID"')
    out, err = capsys.readouterr()
    assert out.strip() == ''

    microhapdb.retrieve.query_mode('population', 'ID == "MHDBP000013"')
    out, err = capsys.readouterr()
    assert '12  MHDBP000013  Jews, Ethiopian  ALFRED' in out


def test_id_mode(capsys):
    microhapdb.retrieve.id_mode(None)
    out, err = capsys.readouterr()
    assert out.strip() == ''

    microhapdb.retrieve.id_mode('NotARealID')
    out, err = capsys.readouterr()
    assert out.strip() == ''

    microhapdb.retrieve.id_mode('rs1363241798')
    out, err = capsys.readouterr()
    line = '21999  MHDBV000022000    GRCh38  chr20  63539694    G,GC  dbSNP151'
    assert line in out


def test_region_mode(capsys):
    with pytest.raises(ValueError) as ve:
        microhapdb.retrieve.region_mode('chr5:1000000-2000000', 'population')
    assert 'region query not supported for table "population"' in str(ve)

    with pytest.raises(ValueError) as ve:
        microhapdb.retrieve.region_mode('chr7:123-456-789', 'locus')
    assert 'cannot parse region "chr7:123-456-789"' in str(ve)

    microhapdb.retrieve.region_mode('chrX', 'variant')
    out, err = capsys.readouterr()
    assert out.strip() == ''

    microhapdb.retrieve.region_mode('chrY', 'locus')
    out, err = capsys.readouterr()
    assert out.strip() == ''

    microhapdb.retrieve.region_mode('chr12:102866940-102866950', None)
    out, err = capsys.readouterr()
    outlines = out.strip().split('\n')
    assert len(outlines) == 8  # 1 locus + 5 variants + 2 header lines
