# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from collections import namedtuple
from files import smartopen
import os
from re import findall, search


Marker = namedtuple('Marker', 'markerid,label,markername,chrom,start,end,variants')
Variant = namedtuple('Variant', 'dbsnpid,chrom,position')


def alfred_locus_list(instream):
    data = instream.read()
    matches = findall(r'(\S+) \| null \| \d+ \| (\S+)', data)
    for locusid, locusname in matches:
        yield locusid, locusname


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
        yield [*values[1:4], 'LOVD', values[0]]


def combine_loci(alfredstream, lovdstream):
    loci = list()
    next(alfredstream)
    for line in alfredstream:
        values = line.strip().split()
        loci.append(values)
    next(lovdstream)
    for line in lovdstream:
        values = line.strip().split()
        loci.append(values)
    loci.sort(key=lambda l: (l[0], int(l[1]), int(l[2])))
    for n, valuelist in enumerate(loci, 1):
        chrom, start, end, source, *idlist = valuelist
        newid = 'MHDBM{:06d}'.format(n)
        yield [newid, 'GRCh38', chrom, start, end, source]
        for xref in idlist:
            yield [xref, 'marker', newid]


def alfred_locus_scrape(locusid, instream, allelefile, variantfile):
    """Scrape a locus detail page for variant and allele data

    We use the locus detail pages to determine which variants are associated
    with each microhap locus, and which compound alleles are represented at
    each locus.
    """
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
    alfredids = [m[4] for m in matches]
    refralleles = [m[1] for m in matches]
    altalleles = [m[2] for m in matches]
    with smartopen(variantfile, 'w') as outstream:
        id_allele_data = zip(dbsnpids, alfredids, refralleles, altalleles)
        for values in id_allele_data:
            print(locusid, *values, sep='\t', file=outstream)

    matches = findall(r'<tr><td>([ACGTDI](-[ACGTDI])+)</td>', data)
    haplotypes = [m[0] for m in matches]
    varstr = ','.join(dbsnpids)
    with smartopen(allelefile, 'w') as outstream:
        for hap in haplotypes:
            hap = hap.replace('-', ',')
            print(locusid, hap, varstr, sep='\t', file=outstream)
