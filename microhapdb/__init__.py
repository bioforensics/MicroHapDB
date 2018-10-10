#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


from microhapdb.util import data_file
from microhapdb import cli
import pandas
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


allelefreqs = pandas.read_table(data_file('allele.tsv'))
loci = pandas.read_table(data_file('locus.tsv'))
populations = pandas.read_table(data_file('population.tsv'))
variants = pandas.read_table(data_file('variant.tsv'))
