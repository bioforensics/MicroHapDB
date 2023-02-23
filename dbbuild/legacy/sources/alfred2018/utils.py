# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import builtins
from collections import defaultdict
from gzip import open as gzopen
import os
from re import findall, search
import sys


def parse_markers_from_table(filename):
    if not os.path.exists(filename):
        return list()
    with smartopen(filename, 'r') as fh:
        return [label for label, xref in alfred_marker_list(fh)]


def smartopen(filename, mode):
    """Smart file handler

    Determines whether to create a compressed or un-compressed file handle
    based on filename extension.
    """
    if mode not in ('r', 'w'):
        raise ValueError('invalid mode "{}"'.format(mode))
    if filename in ['-', None]:
        filehandle = sys.stdin if mode == 'r' else sys.stdout
        return filehandle
    openfunc = builtins.open
    if filename.endswith('.gz'):
        openfunc = gzopen
        mode += 't'
    return openfunc(filename, mode)


def alfred_marker_list(instream):
    data = instream.read()
    matches = findall(r'(\S+) \| null \| \d+ \| (\S+)', data)
    for label, xref in matches:
        yield label, xref


def alfred_marker_detail_scrape(instream):
    """Scrape a marker detail page for variant info"""
    data = instream.read()
    vmatch = search(r'in this (\d+)-site microhap', data)
    assert vmatch, instream.name
    nvariants = int(vmatch.group(1))
    pattern = (
        r'(rs\d+)\s*\((A|C|G|T|Ins|Del|Indel -)/(A|C|G|T|Ins|Del|Indel -)'
        r'( *SNP)*\)[\s\S]*?\(UID= *(\S+)\)'
    )
    matches = findall(pattern, data)
    if len(matches) != nvariants:
        message = 'mismatch: expected {:d} variants'.format(nvariants)
        message += ', only found {:d}'.format(len(matches))
        raise ValueError(message)
    dbsnpids = [m[0] for m in matches]
    match = search('(mh\d+(KK|CP|NK)-\d+)', data)
    assert match, data
    name = match.group(1)
    return name, dbsnpids


def alfred_pop_data(instream, mapping):
    pops_1kgp = set(mapping.ALFRED)
    popdata = dict()
    for line in instream:
        if line.startswith(('----------', 'SI664', 'popName')):
            continue
        values = line.strip().split('\t')
        popmatch = search(r'^([^\(]+)\((\S+)\)', values[0])
        assert popmatch, values[0]
        popname = popmatch.group(1)
        label = popmatch.group(2)
        if label in pops_1kgp:
            continue
        typed_sample_size = int(values[1])
        if label in popdata:
            assert popname == popdata[label]
        else:
            popdata[label] = popname
            yield label, popname, ''
