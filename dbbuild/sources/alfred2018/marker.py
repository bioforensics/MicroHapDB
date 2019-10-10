# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import builtins
from collections import defaultdict, namedtuple
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


def alfred_marker_coords(vcf, mapping):
    rsidcoords = dict()
    for line in vcf:
        if line.startswith('#'):
            continue
        chrnum, posstr, rsid, *values = line.strip.split()
        chrom = chrnum if chrnum.startswith('chr') else 'chr' + chrnum
        pos = int(posstr) - 1
        rsidcoords[rsid] = (chrom, pos)

    for line in mapping:
        if line.startswith('Name'):
            continue
        name, xref, rsidlist = line.strip().split('\t')
        rsids = rsidlist.split(,)
        offsets = [rsidcoords[rsid][1] for rsid in rsids]
        offsets, rsids = zip(*sorted(offsets, rsids))
        chrom = rsidcoords[rsids[0]][0]
        yield name, chrom, offsets, rsids
