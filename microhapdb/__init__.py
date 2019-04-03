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


idmap = pandas.read_csv(data_file('idmap.tsv'), sep='\t')
frequencies = pandas.read_csv(data_file('allele.tsv'), sep='\t')
loci = pandas.read_csv(data_file('locus.tsv'), sep='\t')
populations = pandas.read_csv(data_file('population.tsv'), sep='\t')
variants = pandas.read_csv(data_file('variant.tsv'), sep='\t')
variantmap = pandas.read_csv(data_file('variantmap.tsv'), sep='\t')

tables = {
    'variant': variants,
    'locus': loci,
    'population': populations,
    'allele': frequencies,
}
