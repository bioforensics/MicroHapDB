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


def alfred_marker_coords(vcf, mapping):
    rsidcoords = dict()
    for line in vcf:
        if line.startswith('#'):
            continue
        chrnum, posstr, rsid, *values = line.strip().split()
        chrom = chrnum if chrnum.startswith('chr') else 'chr' + chrnum
        pos = int(posstr) - 1
        rsidcoords[rsid] = (chrom, pos)

    for line in mapping:
        if line.startswith('Name'):
            continue
        name, xref, rsidlist = line.strip().split('\t')
        rsids = rsidlist.split(',')
        offsets = [rsidcoords[rsid][1] for rsid in rsids]
        assert offsets == sorted(offsets)
        chrom = rsidcoords[rsids[0]][0]
        yield name, xref, chrom, offsets, rsids


def alfred_pop_data(instream):
    popdata = dict()
    for line in instream:
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
            yield label, popname


def alfred_frequencies(instream):
    indels = {
        'SI664579L': {'D': 'T', 'I': 'TA'},
        'SI664597L': {'D': 'T', 'I': 'TG'},
        'SI664640A': {'D': 'A', 'I': 'AATAATT'},
    }


    def cleanup(allelestr, markerid, indels):
        allelestr = allelestr.replace('-', ',')
        if 'D' in allelestr or 'I' in allelestr:
            allelestr = allelestr.replace('D', indels[markerid]['D'])
            allelestr = allelestr.replace('I', indels[markerid]['I'])
        return allelestr

    line = next(instream)
    if line.startswith('Created on'):
        next(instream)
    chunks = instream.read().split('-----------------\n')

    for chunk in chunks:
        lines = iter(chunk.split('\n'))
        header1 = next(lines)
        xref, _, _, markerid = header1.split(' | ')
        if '-' not in markerid:
            markerid = markerid[:6] + '-' + markerid[6:]
        header2 = next(lines)
        alleles = header2.split()[2:]
        alleles = [cleanup(a, xref, indels) for a in alleles]
        for line in lines:
            if line.strip() == '':
                continue
            values = line.split('\t')
            popid = search(r'^([^\(]+)\((\S+)\)', values[0]).group(2)
            freqs = values[2:]
            if len(alleles) != len(freqs):
                message = 'WARNING: allele/freq mismatch '
                message += 'for locus ' + markerid
                message += ' and population ' + popid
                message += '; {nfreq} frequencies vs {nall} alleles'.format(
                    nfreq=len(freqs), nall=len(alleles)
                )
                raise ValueError(message)
            for allele, freq in zip(alleles, freqs):
                yield markerid, popid, allele, '{:.4f}'.format(float(freq))
