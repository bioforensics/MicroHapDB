#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import microhapdb


def query_mode(table, querystr):
    if not table:
        raise ValueError('must specify table to invoke a query')
    if table not in microhapdb.tables:
        raise ValueError('unsupported table "{}"'.format(table))
    table = microhapdb.tables[table]
    result = table.query(querystr)
    print(result.to_string())


def id_mode(idstr):
    m = microhapdb.idmap
    prelim = m[(m.XRef == idstr) | (m.mhdbID == idstr)]
    if len(prelim) == 0:
        return
    mhdbids = list(prelim.mhdbID.unique())
    tables = list(prelim.Table.unique())
    assert len(mhdbids) == 1 and len(tables) == 1
    mhdbid = mhdbids[0]
    table = microhapdb.tables[tables[0]]
    result = table[(table.ID == mhdbid)]
    print(result.to_string())


def region_mode(region, table=None):
    chrom, start, end = None, None, None
    if ':' in region:
        chrom, rng = region.split(':')
        if '-' not in rng:
            raise ValueError('cannot parse region "{}"'.format(region))
        startstr, endstr = rng.split('-')
        start, end = int(startstr), int(endstr)
    else:
        chrom = region
    if table not in (None, 'locus', 'variant'):
        msg = 'region query not supported for table "{}"'.format(table)
        raise ValueError(msg)
    if table in (None, 'locus'):
        query = 'Chrom == "{}"'.format(chrom)
        if start is not None:
            query += ' and (Start < {}'.format(end)
            query += ' and End > {})'.format(start)
        result = microhapdb.loci.query(query)
        print(result.to_string())
    if table in (None, 'variant'):
        query = 'Chrom == "{}"'.format(chrom)
        if start is not None:
            query += ' and {} < Position < {}'.format(start, end)
        result = microhapdb.variants.query(query)
        print(result.to_string())
