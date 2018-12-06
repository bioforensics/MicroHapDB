# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from re import search


def population_summary(freqfile):
    """Scrape population data from allele frequency "table" on ALFRED."""
    popdata = dict()
    for line in freqfile:
        if line.startswith(('----------', 'SI664', 'popName')):
            continue
        values = line.strip().split('\t')
        popmatch = search(r'^([^\(]+)\((\S+)\)', values[0])
        assert popmatch, values[0]
        label = popmatch.group(2)
        popname = popmatch.group(1)
        typed_sample_size = int(values[1])
        if label in popdata:
            assert popname == popdata[label]
        else:
            popdata[label] = popname
    for n, label in enumerate(sorted(popdata), 1):
        mhdbid = 'MHDBP{:06d}'.format(n)
        yield mhdbid, label, popdata[label]
