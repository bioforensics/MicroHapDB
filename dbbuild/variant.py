# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from collections import defaultdict, namedtuple
from sys import stderr


Variant = namedtuple('Variant',
                     'dbsnpid,xref,source,markerlabel,chrom,pos,alleles')


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

    def overlap(interval1, interval2):
        if interval1[0] != interval2[0]:
            return False
        return interval1[1] <= interval2[2] and interval1[2] >= interval2[1]

    def merge(intervals):
        tomerge = list()
        for interval in intervals:
            if len(tomerge) > 0:
                current = (tomerge[0][0], tomerge[0][1], tomerge[-1][2])
                if overlap(interval, current):
                    tomerge.append(interval)
                    continue
                else:
                    yield current
                    tomerge = list()
            tomerge.append(interval)
        current = (tomerge[0][0], tomerge[0][1], tomerge[-1][2])
        yield current

    regions = merge(regions)
    regionstrs = ['{:s}:{:d}-{:d}'.format(c, s, e) for c, s, e in regions]
    cmdargs = ['tabix', dbsnpfile, *regionstrs, '>', outfile]
    cmd = ' '.join(cmdargs)
    return cmd


def retrieve_proximal_variants(dbsnp, alfred, lovd, linkoping):
    variants = list()
    variant_labels = defaultdict(list)
    variant_xrefs = defaultdict(list)

    alfred_loaded = set()
    next(alfred)
    for line in alfred:
        values = line.strip().split()
        locuslabel, dbsnpid, xref = values[:3]
        variant_labels[dbsnpid].append(locuslabel)
        variant_xrefs[dbsnpid].append(xref)
        alfred_loaded.add(dbsnpid)

    lovd_loaded = set()
    next(lovd)
    for line in lovd:
        values = line.strip().split()
        dbsnpids = values[5]
        if dbsnpids == '-':
            v = Variant(None, [None], values[-1], [values[0]], values[1],
                        int(values[2]), values[4])
            variants.append(v)
        else:
            locuslabel = values[0]
            for dbsnpid in dbsnpids.split(','):
                variant_labels[dbsnpid].append(locuslabel)
                lovd_loaded.add(dbsnpid)

    linkoping_loaded = set()
    next(linkoping)
    for line in linkoping:
        values = line.strip().split()
        locuslabel = values[0]
        dbsnpid = values[4]
        variant_labels[dbsnpid].append(locuslabel)
        linkoping_loaded.add(dbsnpid)

    alfred_found = set()
    lovd_found = set()
    linkoping_found = set()
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
        if dbsnpid in alfred_loaded:
            alfred_found.add(dbsnpid)
        if dbsnpid in lovd_loaded:
            lovd_found.add(dbsnpid)
        if dbsnpid in linkoping_loaded:
            linkoping_found.add(dbsnpid)

        av = Variant(
            dbsnpid, variant_xrefs[dbsnpid], 'dbSNP151', variant_labels[dbsnpid], chrom,
            position - 1, dbsnpalstr
        )
        variants.append(av)

    variants.sort(key=lambda v: (v.chrom, v.pos))
    for v in variants:
        yield v

    missing = alfred_loaded - alfred_found
    if len(missing) > 0:
        message = '{} variants not found in dbSNP'.format(len(missing))
        message += ': {}'.format(','.join(missing))
        raise ValueError(message)

    missing = lovd_loaded - lovd_found
    if len(missing) > 0:
        message = '{} variants not found in dbSNP'.format(len(missing))
        message += ': {}'.format(','.join(missing))
        raise ValueError(message)

    missing = linkoping_loaded - linkoping_found
    if len(missing) > 0:
        message = '{} variants not found in dbSNP'.format(len(missing))
        message += ': {}'.format(','.join(missing))
        raise ValueError(message)
