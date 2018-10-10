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
    assert len(microhapdb.allelefreqs.Population.unique()) == 84
    # For 13 loci the number of alleles and the number of annotated allele frequencies is mismatched for all populations in the table published by ALFRED.
    # Another locus is only mismatched for a portion of the populations. See `dbbuild/allele-mismatch.log`.
    assert len(microhapdb.allelefreqs.Locus.unique()) == (148 - 13)


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
    16996  SI664616D  SA000077O  T,A,A      0.107
    """
    af = microhapdb.allelefreqs
    assert af.shape == (54637, 4)
    assert len(af.query('Locus == "SI664638H"').Allele.unique()) == 6
    assert len(af.query('Locus == "SI664638H" & Allele == "G,G,G"')) == 83
    freq = af.query('Locus == "SI664638H" & Allele == "G,G,G" & Population == "SA004250L"').Frequency.values[0]
    assert pytest.approx(freq, 0.202)


def test_loci():
    """Microhaplotype locus data

    >>> microhapdb.loci.query('ID == "SI664616D"')
               ID        Name  Chrom     Start       End                          Variants
    43  SI664616D  mh15KK-104     15  63806357  63806495  rs11631544,rs10152453,rs80047978
    >>> microhapdb.loci.query('Chrom == 15 & Start > 50000000 & End < 100000000')
               ID        Name  Chrom     Start       End                          Variants
    43  SI664616D  mh15KK-104     15  63806357  63806495  rs11631544,rs10152453,rs80047978
    46  SI664613A  mh15KK-066     15  52192752  52192827                  rs1063902,rs4219
    """
    loc = microhapdb.loci
    assert loc.shape == (148, 6)
    assert len(loc.query('Chrom == 19')) == 4


def test_populations():
    """Population data

    >>> microhapdb.populations.query('ID == "SA000936S"')
               ID     Name  NumChrom
    48  SA000936S  Koreans       106
    >>> microhapdb.populations.query('Name.str.contains("Afr")')
               ID               Name  NumChrom
    41  SA000101C  African Americans       168
    59  SA004047P  African Americans       122
    75  SA004242M    Afro-Caribbeans       192
    """
    pop = microhapdb.populations
    assert pop.shape == (84, 3)
    assert pop.query('ID == "SA000028K"').Name.values == ['Karitiana']
    assert list(pop.query('Name.str.contains("Jews")').ID.values) == ['SA000015G', 'SA000016H', 'SA000096P', 'SA000490N']


def test_variants():
    """Microhaplotype variant data

    >>> microhapdb.variants.query('ID == "rs80047978"')
                 ID   AlfredID  Chrom     Start       End AlfredAlleles dbSNPAlleles
    298  rs80047978  SI664352A     15  63806494  63806495           A,G          A,G
    >>> microhapdb.variants.query('Chrom == 15 & Start > 50000000 & End < 100000000')
                 ID   AlfredID  Chrom     Start       End AlfredAlleles dbSNPAlleles
    294   rs1063902  SI056461W     15  52192752  52192753           A,C          A,C
    295      rs4219  SI404942X     15  52192826  52192827           G,T          G,T
    296  rs11631544  SI664351Z     15  63806357  63806358           C,T          C,T
    297  rs10152453  SI663904C     15  63806413  63806414           A,C          A,C
    298  rs80047978  SI664352A     15  63806494  63806495           A,G          A,G
    """
    v = microhapdb.variants
    assert v.shape == (405, 7)
    assert len(v.query('Chrom == 12')) == 15
