# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from collections import namedtuple
from sys import stderr


Variant = namedtuple('Variant', 'dbsnpid,alfredid,label,chrom,pos,alleles')


def retrieve_proximal_variants(alfred, dbsnp):
    alfred_found = set()
    alfred_data = dict()
    next(alfred)
    for line in alfred:
        values = line.strip().split()
        dbsnpid = values[1]
        alfred_data[dbsnpid] = values
    next(dbsnp)
    for line in dbsnp:
        dbsnpvals = line.strip().split()
        dbsnpid = dbsnpvals[2]
        dba = dbsnpvals[3:5]
        dbsnpalleles = set([n for v in dba for n in v.split(',')])
        dbsnpalstr = ','.join(sorted(dbsnpalleles))
        chrom = dbsnpvals[0]
        position = int(dbsnpvals[1])
        alfredid = None
        locuslabel = None
        if dbsnpid in alfred_data:
            alfred_found.add(dbsnpid)
            alfredid = alfred_data[dbsnpid][2]
            locuslabel = alfred_data[dbsnpid][0]
        av = Variant(dbsnpid, alfredid, locuslabel, chrom, position - 1,
                     dbsnpalstr)
        yield av
    missing = set(alfred_data.keys()) - alfred_found
    if len(missing) > 0:
        message = '{} variants not found in dbSNP'.format(len(missing))
        message += ': {}'.format(','.join(missing))
        raise ValueError(message)
