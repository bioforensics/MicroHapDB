# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from collections import defaultdict
import microhapdb


def standardize_ids(idents):
    result = microhapdb.populations[
        (microhapdb.populations.ID.isin(idents)) | (microhapdb.populations.Name.isin(idents))
    ]
    return result.ID


def print_table(table, **kwargs):
    print(table.to_string(index=False))


def print_detail(table):
    for n, row in table.iterrows():
        print('----------------------------------------------------------[ MicroHapulator ]----')
        print('{name:s}    ({id:s}; source={src:s})\n'.format(id=row.ID, name=row.Name, src=row.Source))
        result = microhapdb.frequencies[microhapdb.frequencies.Population == row.ID]
        markers = sorted(result.Marker.unique())
        frequencies = len(result.Frequency)
        print('- {freq:d} total allele frequencies available\n  for {mark:d} markers'.format(
            mark=len(markers), freq=frequencies
        ))
        counter = defaultdict(int)
        for marker in markers:
            alleles = result[result.Marker == marker]
            nalleles = len(alleles.Allele)
            counter[nalleles] += 1
        print('\n# Alleles | # Markers\n---------------------')
        for nalleles, nmarkers in sorted(counter.items(), reverse=True):
            print('       {na:3d}|{nm:s}'.format(na=nalleles, nm='*' * nmarkers))
        print('--------------------------------------------------------------------------------\n')
