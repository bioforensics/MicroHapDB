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


def load_variants():
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
    return pandas.read_table(microhapdb.data_file('variant.tsv'))


def test_variants():
    v = microhapdb.variants
    assert v.shape == (405, 7)
    assert len(v.query('Chrom == 12')) == 15


variants = load_variants()
