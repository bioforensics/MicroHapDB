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


def allele_positions(locusid):
    """Convenience function for grabbing a locus' genomic coordinates.

    Loci can be accessed by their ALFRED IDs (SI...) or their names (mh...).
    """
    q = 'ID == "{id}" | Name == "{id}"'.format(id=locusid)
    locus = loci.query(q)
    assert len(locus) == 1
    dbsnpids = locus['Variants'].iloc[0].split(',')
    var = variants[variants['ID'].isin(dbsnpids)]
    for index, row in var.iterrows():
        yield row['Chrom'], row['Start'], row['End']
