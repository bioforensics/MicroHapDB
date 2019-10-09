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


Marker = namedtuple('Marker', 'markerid,label,markername,chrom,start,end,variants')
Variant = namedtuple('Variant', 'dbsnpid,chrom,position')


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
    return dbsnpids


def alfred_loci_prelim(freqstream, coordstream):
    locusids = dict()
    for locusid, locusname in alfred_locus_list(freqstream):
        locusids[locusname] = locusid
    next(coordstream)
    for line in coordstream:
        values = line.strip().split()
        locusid = locusids[values[0]]
        yield [*values[1:], 'ALFRED', locusid, values[0]]


def lovd_loci_prelim(instream):
    for line in instream:
        if line.startswith(' '):
            continue
        values = line.strip().split(',')
        offsets = [int(o) for o in values[-1].split(':')]
        start = int(values[2])
        values[2] = start + offsets[0]
        values[3] = start + offsets[-1] + 1
        yield [*values[1:4], 'doi:10.1016/j.fsigen.2018.05.008', values[0]]


def linköping_loci_prelim(defstream, varstream):
    coords = dict()
    for line in varstream:
        if line.startswith('#'):
            continue
        chrnum, pos, rsid, *data = line.split('\t')
        chrom = 'chr' + chrnum
        coords[rsid] = (chrnum, int(pos) - 1)

    microhaps = defaultdict(list)
    for line in defstream:
        if not line.startswith('rs'):
            continue
        rsid, label = line.strip().split('\t')
        microhaps[label].append(rsid)

    for label, snplist in microhaps.items():
        loclist = [coords[rsid] for rsid in snplist]
        chrom = set([loc[0] for loc in loclist])
        assert len(chrom) == 1, chrom
        chromstr = 'chr{:s}'.format(list(chrom)[0])
        coordlist = [loc[1] for loc in loclist]
        start = min(coordlist)
        end = max(coordlist) + 1
        yield chromstr, start, end, 'ISFG2019:P597', label


def combine_loci(alfredstream, lovdstream, linkstream):
    loci = list()
    for varstream in (alfredstream, lovdstream, linkstream):
        next(varstream)
        for line in varstream:
            values = line.strip().split()
            loci.append(values)
    loci.sort(key=lambda l: (l[0], int(l[1]), int(l[2])))
    for n, valuelist in enumerate(loci, 1):
        chrom, start, end, source, *idlist = valuelist
        newid = 'MHDBM{:06d}'.format(n)
        yield [newid, 'GRCh38', chrom, start, end, source]
        for xref in idlist:
            yield [xref, 'marker', newid]
