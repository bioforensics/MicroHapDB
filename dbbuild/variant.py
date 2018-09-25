# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from collections import namedtuple


Variant = namedtuple('Variant', 'dbsnpid, chrom, position, refr, alt')
AggVariant = namedtuple(
    'AggVariant',
    'dbsnpid, alfredid, chrom, start, end, alfalleles, dbsnpalleles'
)


def variant_xref(alfred, dbsnp):
    dbsnpids = set()
    next(alfred)
    for line in alfred:
        dbsnpid = line.split()[1]
        dbsnpids.add(dbsnpid)
    dbsnpids_found = set()
    for line in dbsnp:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        dbsnpid = fields[2]
        if dbsnpid not in dbsnpids:
            continue
        chrom = fields[0]
        position = int(fields[1])
        refr = fields[3]
        alt = fields[4]
        yield Variant(dbsnpid, chrom, position, refr, alt)
        dbsnpids_found.add(dbsnpid)
    if dbsnpids > dbsnpids_found:
        missing = sorted(dbsnpids - dbsnpids_found)
        message = '{} dbSNP IDs missing: {}'.format(len(missing), missing)
        raise ValueError(message)


def combine_variants(alfred, dbsnp):
    alfred_data = dict()
    next(alfred)
    for line in alfred:
        values = line.strip().split()
        dbsnpid = values[1]
        alfred_data[dbsnpid] = values
    next(dbsnp)
    for line in dbsnp:
        dbsnpvals = line.strip().split()
        dbsnpid = dbsnpvals[0]
        alfredvals = alfred_data[dbsnpid]
        alfredid = alfredvals[2]
        chrom = dbsnpvals[1]
        alfalleles = set(alfredvals[3:5])
        alfalstr = ','.join(sorted(alfalleles))
        dba = dbsnpvals[3:5]
        dbsnpalleles = set([n for v in dba for n in v.split(',')])
        dbsnpalstr = ','.join(sorted(dbsnpalleles))
        position = int(dbsnpvals[2])
        av = AggVariant(dbsnpid, alfredid, chrom, position - 1, position,
                        alfalstr, dbsnpalstr)
        yield av
