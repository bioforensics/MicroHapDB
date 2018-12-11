# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from collections import namedtuple
from sys import stderr


Variant = namedtuple('Variant',
                     'dbsnpid,xref,source,locuslabel,chrom,pos,alleles')


def dbsnp_subset_command(locusfile, dbsnpfile, outfile):
    regions = list()
    with open(locusfile, 'r') as instream:
        next(instream)
        for line in instream:
            values = line.strip().split()
            chrom = values[2][3:]
            start, end = int(values[3]), int(values[4])
            length = start - end
            diff = 500 - length
            extend = diff / 2
            newstart = int(start - extend)
            newend = int(end + extend)
            regions.append((chrom, newstart, newend))
    regions = sorted(regions, key=lambda r: (r[0], r[1], r[2]))
    regionstrs = ['{:s}:{:d}-{:d}'.format(c, s, e) for c, s, e in regions]
    cmdargs = ['tabix', dbsnpfile, *regionstrs, '>', outfile]
    cmd = ' '.join(cmdargs)
    return cmd


def retrieve_proximal_variants(dbsnp, alfred, lovd):
    variants = list()

    alfred_found = set()
    alfred_data = dict()
    next(alfred)
    for line in alfred:
        values = line.strip().split()
        dbsnpid = values[1]
        alfred_data[dbsnpid] = values

    lovd_found = set()
    lovd_data = dict()
    next(lovd)
    for line in lovd:
        values = line.strip().split()
        dbsnpids = values[5]
        if dbsnpids == '-':
            v = Variant(None, None, values[-1], values[0], values[1],
                        int(values[2]) - 1, values[4])
            variants.append(v)
        else:
            for dbsnpid in dbsnpids.split(','):
                lovd_data[dbsnpid] = values

    for line in dbsnp:
        dbsnpvals = line.strip().split()
        dbsnpid = dbsnpvals[2]
        dba = dbsnpvals[3:5]
        dbsnpalleles = set([n for v in dba for n in v.split(',')])
        dbsnpalstr = ','.join(sorted(dbsnpalleles))
        chrom = dbsnpvals[0]
        if not chrom.startswith('chr'):
            chrom = 'chr' + chrom
        position = int(dbsnpvals[1])
        xref = None
        locuslabel = None
        if dbsnpid in alfred_data:
            alfred_found.add(dbsnpid)
            xref = alfred_data[dbsnpid][2]
            locuslabel = alfred_data[dbsnpid][0]
        elif dbsnpid in lovd_data:
            lovd_found.add(dbsnpid)
            locuslabel = lovd_data[dbsnpid][0]
        av = Variant(dbsnpid, xref, 'dbSNP151', locuslabel, chrom, position - 1,
                     dbsnpalstr)
        variants.append(av)

    variants.sort(key=lambda v: (v.chrom, v.pos))
    for v in variants:
        yield v

    missing = set(alfred_data.keys()) - alfred_found
    if len(missing) > 0:
        message = '{} variants not found in dbSNP'.format(len(missing))
        message += ': {}'.format(','.join(missing))
        raise ValueError(message)
