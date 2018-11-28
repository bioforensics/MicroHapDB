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
    assert len(microhapdb.allelefreqs.Population.unique()) == 96
    # For many loci there is a discrepancy between the number of annotated alleles and allele frequencies in the table published by ALFRED.
    # See `dbbuild/allele-mismatch.log` for details.
    # The problematic allele frequencies are discarded, and only the congruent data are retained.
    # For 1 locus, however, there is a discrepancy for every single single data point, and thus no frequency data for that locus is retained.
    assert len(microhapdb.allelefreqs.Locus.unique()) == (198 - 1)


def test_allele_frequencies():
    """Allele frequencies for 148 loci across 84 populations

    >>> alleles = microhapdb.allelefreqs.query('Locus == "SI664616D"').Allele
    >>> list(alleles.unique())
    ['T,A,G', 'T,A,A', 'T,C,G', 'C,A,G']
    >>> pop_freqs = microhapdb.allelefreqs.query('Locus == "SI664616D" & Allele == "T,A,A"')
    >>> len(pop_freqs)
    83
    >>> microhapdb.allelefreqs.query('Locus == "SI664616D" & Allele == "T,A,A" & Population == "SA000077O"')
               Locus Population Allele  Frequency
    23915  SI664616D  SA000077O  T,A,A      0.107
    """
    af = microhapdb.allelefreqs
    assert af.shape == (71535, 4)
    assert len(af.query('Locus == "SI664638H"').Allele.unique()) == 6
    assert len(af.query('Locus == "SI664638H" & Allele == "G,G,G"')) == 83
    freq = af.query('Locus == "SI664638H" & Allele == "G,G,G" & Population == "SA004250L"').Frequency.values[0]
    assert pytest.approx(freq, 0.202)


def test_loci():
    """Microhaplotype locus data

    >>> microhapdb.loci.query('ID == "SI664616D"')
               ID        Name  Chrom     Start       End                          Variants
    76  SI664616D  mh15KK-104     15  63806357  63806495  rs11631544,rs10152453,rs80047978
    >>> microhapdb.loci.query('Chrom == 15 & Start > 50000000 & End < 100000000')
               ID        Name  Chrom     Start       End                          Variants
    73  SI664613A  mh15KK-066     15  52192752  52192827                  rs1063902,rs4219
    76  SI664616D  mh15KK-104     15  63806357  63806495  rs11631544,rs10152453,rs80047978
    """
    loc = microhapdb.loci
    assert loc.shape == (198, 6)
    assert len(loc.query('Chrom == 19')) == 5


def test_populations():
    """Population data

    >>> microhapdb.populations.query('ID == "SA000936S"')
               ID     Name  NumChrom
    47  SA000936S  Koreans       106
    >>> microhapdb.populations.query('Name.str.contains("Afr")')
               ID               Name  NumChrom
    40  SA000101C  African Americans       168
    58  SA004047P  African Americans       122
    74  SA004242M    Afro-Caribbeans       192
    """
    pop = microhapdb.populations
    assert pop.shape == (83, 3)
    assert pop.query('ID == "SA000028K"').Name.values == ['Karitiana']
    assert list(pop.query('Name.str.contains("Jews")').ID.values) == ['SA000015G', 'SA000016H', 'SA000096P', 'SA000490N']


def test_variants():
    """Microhaplotype variant data

    >>> microhapdb.variants.query('ID == "rs80047978"')
                 ID   AlfredID  Chrom     Start       End AlfredAlleles dbSNPAlleles
    422  rs80047978  SI664352A     15  63806494  63806495           A,G          A,G
    >>> microhapdb.variants.query('Chrom == 15 & Start > 50000000 & End < 100000000')
                 ID   AlfredID  Chrom     Start       End AlfredAlleles dbSNPAlleles
    418   rs1063902  SI056461W     15  52192752  52192753           A,C          A,C
    419      rs4219  SI404942X     15  52192826  52192827           G,T          G,T
    420  rs11631544  SI664351Z     15  63806357  63806358           C,T          C,T
    421  rs10152453  SI663904C     15  63806413  63806414           A,C          A,C
    422  rs80047978  SI664352A     15  63806494  63806495           A,G          A,G
    """
    v = microhapdb.variants
    assert v.shape == (559, 7)
    assert len(v.query('Chrom == 12')) == 20


def test_allele_positions():
    """
    >>> list(microhapdb.allele_positions('mh19KK-301'))
    [(19, 50938487, 50938488), (19, 50938502, 50938503), (19, 50938526, 50938527), (19, 50938550, 50938551)]
    >>> list(microhapdb.allele_positions('SI664193D'))
    [(19, 4852124, 4852125), (19, 4852324, 4852325)]
    """
    pass
