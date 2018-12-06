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
from microhapdb import retrieve
from microhapdb.retrieve import id_xref
import pandas
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


idmap = pandas.read_table(data_file('idmap.tsv'))
frequencies = pandas.read_table(data_file('allele.tsv'))
loci = pandas.read_table(data_file('locus.tsv'))
populations = pandas.read_table(data_file('population.tsv'))
variants = pandas.read_table(data_file('variant.tsv'))
variantmap = pandas.read_table(data_file('variantmap.tsv'))

tables = {
    'variant': variants,
    'locus': loci,
    'population': populations,
    'allele': frequencies,
}
