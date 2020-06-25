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
    columns = ['Name', 'PermID', 'Reference', 'Chrom', 'Offsets', 'Ae', 'In', 'Fst', 'Source']
    if popid is None:
        defaults = pandas.read_csv(data_file('marker.tsv'), sep='\t')
        defaults = defaults[['Name', 'Ae']]
        markers = markers.drop(columns=['Ae']).join(defaults.set_index('Name'), on='Name')[columns]
    else:
        aes = pandas.read_csv(data_file('marker-aes.tsv'), sep='\t')
        if popid not in aes.Population.unique():
            raise ValueError(f'no Ae data for population "{popid}"')
        popaes = aes[aes.Population == popid].drop(columns=['Population'])
        markers = markers.drop(columns=['Ae']).join(popaes.set_index('Marker'), on='Name')[columns]


def set_reference(refr):
    global markers
    assert refr in (37, 38)
    columns = ['Name', 'PermID', 'Reference', 'Chrom', 'Offsets', 'Ae', 'In', 'Fst', 'Source']
    if refr == 38:
        defaults = pandas.read_csv(data_file('marker.tsv'), sep='\t')[['Name', 'Reference', 'Offsets']]
        markers = markers.drop(columns=['Reference', 'Offsets']).join(defaults.set_index('Name'), on='Name')[columns]
    else:
        o37 = pandas.read_csv(data_file('marker-offsets-GRCh37.tsv'), sep='\t')
        markers = markers.drop(columns=['Reference', 'Offsets']).join(o37.set_index('Marker'), on='Name')[columns]


markers = pandas.read_csv(data_file('marker.tsv'), sep='\t')
populations = pandas.read_csv(data_file('population.tsv'), sep='\t')
frequencies = pandas.read_csv(data_file('frequency.tsv'), sep='\t')
variantmap = pandas.read_csv(data_file('variantmap.tsv'), sep='\t')
idmap = pandas.read_csv(data_file('idmap.tsv'), sep='\t')
sequences = pandas.read_csv(data_file('sequences.tsv'), sep='\t')
indels = pandas.read_csv(data_file('indels.tsv'), sep='\t')
