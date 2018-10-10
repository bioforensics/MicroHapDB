#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


import pandas
import pytest


@pytest.fixture(autouse=True, scope='session')
def pandas_terminal_width():
    pandas.set_option('display.width', 1000)
    pandas.set_option('display.max_columns', 1000)
