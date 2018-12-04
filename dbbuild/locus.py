# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from collections import namedtuple
from files import sopen, tmpfile
import os
from re import findall, search


Locus = namedtuple('Locus', 'locusid,label,locusname,chrom,start,end,variants')


# Populate a list of all microhap loci, but only if the locus table has already
# been created from the allele frequency database.
ALL_LOCI = list()
if os.path.isfile(tmpfile('microhap-loci-prelim.tsv')):
    with open(tmpfile('microhap-loci-prelim.tsv')) as instream:
        next(instream)
        for line in instream:
            locusid, name = line.strip().split()
            ALL_LOCI.append(locusid)


def locus_list(instream):
    data = instream.read()
    matches = findall(r'(\S+) \| null \| \d+ \| (\S+)', data)
    for locusid, locusname in matches:
        yield locusid, locusname


def locus_scrape(locusid, instream, allelefile, variantfile):
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
    with sopen(variantfile, 'w') as outstream:
        id_allele_data = zip(dbsnpids, alfredids, refralleles, altalleles)
        for values in id_allele_data:
            print(locusid, *values, sep='\t', file=outstream)

    matches = findall(r'<tr><td>([ACGTDI](-[ACGTDI])+)</td>', data)
    haplotypes = [m[0] for m in matches]
    varstr = ','.join(dbsnpids)
    with sopen(allelefile, 'w') as outstream:
        for hap in haplotypes:
            hap = hap.replace('-', ',')
            print(locusid, hap, varstr, sep='\t', file=outstream)


def combine_locus_data(idstream, allelestream, variantstream):
    """Aggregate locus data from ALFRED and dbSNP."""
    variants = dict()
    next(variantstream)
    for line in variantstream:
        values = line.strip().split()
        vid = values[0]
        variants[vid] = values

    locusvariants = dict()
    next(allelestream)
    for line in allelestream:
        values = line.strip().split()
        locusid = values[0]
        variantids = values[2]
        locusvariants[locusid] = variantids

    next(idstream)
    for n, line in enumerate(idstream, 1):
        locusid = 'MHDBL{:06d}'.format(n)
        label, locusname = line.strip().split()
        varstr = locusvariants[label]
        varlist = varstr.split(',')
        vardata = [variants[v] for v in varlist]
        chrom = vardata[0][2]
        start = min([int(d[3]) for d in vardata])
        end = max([int(d[4]) for d in vardata])
        yield Locus(locusid, label, locusname, chrom, start, end, varlist)
