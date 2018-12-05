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
    >>> freqs.query('Locus == "MHDBL000077" and Allele == "T,A,A" '
    ...             'and Population == "MHDBP000034"')
                 Locus   Population Allele  Frequency
    30638  MHDBL000077  MHDBP000034  T,A,A      0.107
    """
    af = microhapdb.frequencies
    assert af.shape == (82167, 4)
    assert len(af.query('Locus == "MHDBL000135"').Allele.unique()) == 6
    assert len(af.query('Locus == "MHDBL000135" and Allele == "G,G,G"')) == 83
    query = ('Locus == "MHDBL000135" and Allele == "G,G,G" '
             'and Population == "MHDBP000092"')
    freq = af.query(query).Frequency.values[0]
    assert pytest.approx(freq, 0.202)


def test_loci():
    """Microhaplotype locus data

    >>> microhapdb.loci.query('ID == "MHDBL000077"')
                 ID Reference  Chrom     Start       End  Source
    76  MHDBL000077    GRCh38  chr15  63806358  63806495  ALFRED
    """
    loc = microhapdb.loci
    assert loc.shape == (198, 6)
    assert len(loc.query('Chrom == "chr19"')) == 5
#
# FIXME NEED MORE LOGIC AND LESS QUERIES!
#
# def test_populations():
#     """Population data
#
#     >>> microhapdb.populations.query('ID == "SA000936S"')
#                ID     Name  NumChrom
#     47  SA000936S  Koreans       106
#     >>> microhapdb.populations.query('Name.str.contains("Afr")')
#                ID               Name  NumChrom
#     40  SA000101C  African Americans       168
#     58  SA004047P  African Americans       122
#     74  SA004242M    Afro-Caribbeans       192
#     """
#     pop = microhapdb.populations
#     assert pop.shape == (83, 3)
#     assert pop.query('ID == "SA000028K"').Name.values == ['Karitiana']
#     assert list(pop.query('Name.str.contains("Jews")').ID.values) == ['SA000015G', 'SA000016H', 'SA000096P', 'SA000490N']
#
#
# def test_variants():
#     """Microhaplotype variant data
#
#     >>> microhapdb.variants.query('ID == "rs80047978"')
#                  ID   AlfredID  Chrom     Start       End AlfredAlleles dbSNPAlleles
#     422  rs80047978  SI664352A     15  63806494  63806495           A,G          A,G
#     >>> microhapdb.variants.query('Chrom == 15 & Start > 50000000 & End < 100000000')
#                  ID   AlfredID  Chrom     Start       End AlfredAlleles dbSNPAlleles
#     418   rs1063902  SI056461W     15  52192752  52192753           A,C          A,C
#     419      rs4219  SI404942X     15  52192826  52192827           G,T          G,T
#     420  rs11631544  SI664351Z     15  63806357  63806358           C,T          C,T
#     421  rs10152453  SI663904C     15  63806413  63806414           A,C          A,C
#     422  rs80047978  SI664352A     15  63806494  63806495           A,G          A,G
#     """
#     v = microhapdb.variants
#     assert v.shape == (559, 7)
#     assert len(v.query('Chrom == 12')) == 20
#
#
# def test_allele_positions():
#     """
#     >>> list(microhapdb.allele_positions('mh19KK-301'))
#     [(19, 50938487, 50938488), (19, 50938502, 50938503), (19, 50938526, 50938527), (19, 50938550, 50938551)]
#     >>> list(microhapdb.allele_positions('SI664193D'))
#     [(19, 4852124, 4852125), (19, 4852324, 4852325)]
#     """
#     pass
#
