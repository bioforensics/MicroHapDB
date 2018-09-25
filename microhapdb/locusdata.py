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


def load_loci():
    """Microhaplotype locus data

    >>> microhapdb.loci.query('ID == "SI664616D"')
               ID        Name  Chrom     Start       End                          Variants
    43  SI664616D  mh15KK-104     15  63806357  63806495  rs11631544,rs10152453,rs80047978
    >>> microhapdb.loci.query('Chrom == 15 & Start > 50000000 & End < 100000000')
               ID        Name  Chrom     Start       End                          Variants
    43  SI664616D  mh15KK-104     15  63806357  63806495  rs11631544,rs10152453,rs80047978
    46  SI664613A  mh15KK-066     15  52192752  52192827                  rs1063902,rs4219
    """
    return pandas.read_table(microhapdb.data_file('locus.tsv'))


def test_loci():
    loc = microhapdb.loci
    assert loc.shape == (148, 6)
    assert len(loc.query('Chrom == 19')) == 4


loci = load_loci()
