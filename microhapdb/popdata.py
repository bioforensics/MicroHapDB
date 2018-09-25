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


def load_populations():
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
    return pandas.read_table(microhapdb.data_file('population.tsv'))


def test_populations():
    pop = microhapdb.populations
    assert pop.shape == (84, 3)
    assert pop.query('ID == "SA000028K"').Name.values == ['Karitiana']
    assert list(pop.query('Name.str.contains("Jews")').ID.values) == ['SA000015G', 'SA000016H', 'SA000096P', 'SA000490N']


populations = load_populations()
