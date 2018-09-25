#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import microhapdb
import pandas
import pytest


def load_allele_frequencies():
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
    return pandas.read_table(microhapdb.data_file('allele.tsv'))


def test_allele_frequencies():
    af = microhapdb.allelefreqs
    assert af.shape == (54637, 4)
    assert len(af.query('Locus == "SI664638H"').Allele.unique()) == 6
    assert len(af.query('Locus == "SI664638H" & Allele == "G,G,G"')) == 83
    freq = af.query('Locus == "SI664638H" & Allele == "G,G,G" & Population == "SA004250L"').Frequency.values[0]
    assert pytest.approx(freq, 0.202)


allelefreqs = load_allele_frequencies()
