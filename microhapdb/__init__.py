#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


from microhapdb.util import data_file
from microhapdb.alfreqdata import allelefreqs
from microhapdb.locusdata import loci
from microhapdb.popdata import populations
from microhapdb.variantdata import variants
from microhapdb import cli


from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
