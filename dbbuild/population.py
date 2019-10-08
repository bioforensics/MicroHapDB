# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from re import search


def population_summary(alfredstream, lovdstream):
    popdata = dict()
    for line in alfredstream:
        if line.startswith(('----------', 'SI664', 'popName')):
            continue
        values = line.strip().split('\t')
        popmatch = search(r'^([^\(]+)\((\S+)\)', values[0])
        assert popmatch, values[0]
        label = popmatch.group(2)
        popname = popmatch.group(1)
        typed_sample_size = int(values[1])
        if label in popdata:
            assert popname == popdata[label][0]
        else:
            popdata[label] = (popname, 'ALFRED')
    for line in lovdstream:
        label, source = line.strip().split()
        popdata[label] = (label, 'LOVD')
    pops = list(popdata.items())
    pops.append(('Swedish', ('Swedish', 'Linköping')))
    pops = sorted(pops, key=lambda d: d[1][0])
    for n, (label, (name, source)) in enumerate(pops, 1):
        mhdbid = 'MHDBP{:06d}'.format(n)
        yield mhdbid, label, name, source
