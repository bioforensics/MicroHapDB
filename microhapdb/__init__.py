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
from microhapdb import marker
from microhapdb import panel
from microhapdb import population
import os
from pkg_resources import resource_filename
import pandas
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


def set_ae_population(popid=None):
    global markers
    markers = pandas.read_csv(data_file('marker.tsv'), sep='\t')
    if popid is None:
        return
    else:
        aes = pandas.read_csv(data_file('marker-aes.tsv'), sep='\t')
        if popid not in aes.Population.unique():
            raise ValueError(f'no Ae data for population "{popid}"')
        popaes = aes[aes.Population == popid].drop(columns=['Population'])
        columns = ['Name', 'PermID', 'Reference', 'Chrom', 'Offsets', 'Ae', 'In', 'Fst', 'Source']
        markers = markers.drop(columns=['Ae']).join(popaes.set_index('Marker'), on='Name')[columns]


def set_reference(refr):
    global markers
    assert refr in (37, 38)
    markers = pandas.read_csv(data_file('marker.tsv'), sep='\t')
    if refr == 38:
        return
    o37 = pandas.read_csv(data_file('marker-offsets-GRCh37.tsv'), sep='\t')
    columns = ['Name', 'PermID', 'Reference', 'Chrom', 'Offsets', 'Ae', 'In', 'Fst', 'Source']
    markers = markers.drop(columns=['Reference', 'Offsets']).join(o37.set_index('Marker'), on='Name')[columns]



markers = pandas.read_csv(data_file('marker.tsv'), sep='\t')
populations = pandas.read_csv(data_file('population.tsv'), sep='\t')
frequencies = pandas.read_csv(data_file('frequency.tsv'), sep='\t')
variantmap = pandas.read_csv(data_file('variantmap.tsv'), sep='\t')
idmap = pandas.read_csv(data_file('idmap.tsv'), sep='\t')
sequences = pandas.read_csv(data_file('sequences.tsv'), sep='\t')
indels = pandas.read_csv(data_file('indels.tsv'), sep='\t')
