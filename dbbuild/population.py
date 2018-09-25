# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from collections import namedtuple
from re import search


Population = namedtuple('Population', 'popid, popname, numchrom')


def population_summary(freqfile):
    """Scrape population data from allele frequency "table" on ALFRED."""
    popdata = list()
    for line in freqfile:
        if not line.startswith('popName'):
            continue
        for line in freqfile:
            if line.startswith('--------'):
                break
            values = line.strip().split('\t')
            popmatch = search(r'^([^\(]+)\((\S+)\)', values[0])
            assert popmatch, values[0]
            popid = popmatch.group(2)
            popname = popmatch.group(1)
            numchrom = int(values[1])
            pop = Population(popid, popname, numchrom)
            popdata.append(pop)
        break
    for pop in sorted(popdata, key=lambda p: p.popid):
        yield pop
