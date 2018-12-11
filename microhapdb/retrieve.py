#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import microhapdb


def fetch_by_query(table, querystr):
    """Retrieve data using a Pandas-style query.

    >>> for result in fetch_by_query('locus', 'Chrom == "chr8"'): print(result)
                  ID Reference Chrom     Start       End   AvgAe  Source
    202  MHDBL000203    GRCh38  chr8   1194331   1194395  2.6273    LOVD
    203  MHDBL000204    GRCh38  chr8   3659270   3659482  3.9708  ALFRED
    204  MHDBL000205    GRCh38  chr8  11738320  11738460  2.2865  ALFRED
    >>> querystr = 'Locus == "MHDBL000177" and Population == "MHDBP000092"'
    >>> for result in fetch_by_query('allele', querystr): print(result)
                 Locus   Population Allele  Frequency
    66015  MHDBL000177  MHDBP000092    C,C      0.661
    66016  MHDBL000177  MHDBP000092    C,T      0.339
    66017  MHDBL000177  MHDBP000092    T,C      0.000
    66018  MHDBL000177  MHDBP000092    T,T      0.000
    """
    if not table:
        raise ValueError('must specify table to invoke a query')
    if table not in microhapdb.tables:
        raise ValueError('unsupported table "{}"'.format(table))
    table = microhapdb.tables[table]
    result = table.query(querystr)
    if len(result) > 0:
        yield result


def fetch_by_id(idvalue):
    """Retrieve data using any internal/external ID, name, or label.

    The main data tables are indexed using internal MicroHapDB IDs, but we
    maintain an auxiliary table mapping names, IDs, and labels from source
    databases to the internal IDs. Using this table, we can retrieve the
    relevant record(s) using any valid record identifier.

    >>> for result in fetch_by_id('SI605775E'): print(result)
                       ID Reference  Chrom  Position Alleles    Source
    10482  MHDBV000010483    GRCh38  chr13  50313423     C,T  dbSNP151
    >>> for result in fetch_by_id('rs690302'): print(result)
                       ID Reference  Chrom  Position  Alleles    Source
    19301  MHDBV000019302    GRCh38  chr18   8892896  A,C,G,T  dbSNP151
    >>> for result in fetch_by_id('mh06PK-25713'): print(result)
                  ID Reference Chrom     Start       End   AvgAe Source
    192  MHDBL000193    GRCh38  chr6  31196947  31197015  3.0005   LOVD
    >>> for result in fetch_by_id('mh19CP-007'): print(result)
                  ID Reference  Chrom     Start       End   AvgAe  Source
    109  MHDBL000110    GRCh38  chr19  14310740  14310781  3.0813  ALFRED
    >>> for result in fetch_by_id('SA004109O'): print(result)
                 ID       Name  Source
    15  MHDBP000016  Colombian  ALFRED
    """
    m = microhapdb.idmap
    result = m[(m.XRef == idvalue) | (m.mhdbID == idvalue)]
    if len(result) == 0:
        return
    tables = list(result.Table.unique())
    assert len(tables) == 1, tables
    table = microhapdb.tables[str(tables[0])]
    yield table[table.ID.isin(result.mhdbID)]


def fetch_by_region(region, table=None):
    """Retrieve data using a genomic range.

    Optionally, use `table` to restrict the results to only variants or only
    loci.

    >>> for result in fetch_by_region('chr13', table='locus'): print(result)
                 ID Reference  Chrom      Start        End   AvgAe  Source
    55  MHDBL000056    GRCh38  chr13   23191402   23191542  3.6440  ALFRED
    56  MHDBL000057    GRCh38  chr13   24343963   24343994  3.0655  ALFRED
    57  MHDBL000058    GRCh38  chr13   46291795   46291987  4.0035  ALFRED
    58  MHDBL000059    GRCh38  chr13   50313424   50313589  2.3932  ALFRED
    59  MHDBL000060    GRCh38  chr13   53486692   53486837  6.0444  ALFRED
    60  MHDBL000061    GRCh38  chr13   66138600   66138696  3.4361  ALFRED
    61  MHDBL000062    GRCh38  chr13   94894396   94894513  2.0150  ALFRED
    62  MHDBL000063    GRCh38  chr13  110154352  110154505  3.7978  ALFRED
    >>> for result in fetch_by_region('chr18:24557400-24557450'): print(result)
                  ID Reference  Chrom     Start       End  AvgAe  Source
    105  MHDBL000106    GRCh38  chr18  24557355  24557490  2.658  ALFRED
                       ID Reference  Chrom  Position Alleles    Source
    19445  MHDBV000019446    GRCh38  chr18  24557400     C,T  dbSNP151
    19446  MHDBV000019447    GRCh38  chr18  24557402     A,C  dbSNP151
    19447  MHDBV000019448    GRCh38  chr18  24557414     A,T  dbSNP151
    19448  MHDBV000019449    GRCh38  chr18  24557416     A,G  dbSNP151
    19449  MHDBV000019450    GRCh38  chr18  24557431     A,G  dbSNP151
    19450  MHDBV000019451    GRCh38  chr18  24557443     A,G  dbSNP151
    19451  MHDBV000019452    GRCh38  chr18  24557447   C,G,T  dbSNP151
    19452  MHDBV000019453    GRCh38  chr18  24557448     A,G  dbSNP151
    """
    chrom, start, end = None, None, None
    if ':' in region:
        chrom, rng = region.split(':')
        if rng.count('-') != 1:
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
        if len(result) > 0:
            yield result
    if table in (None, 'variant'):
        query = 'Chrom == "{}"'.format(chrom)
        if start is not None:
            query += ' and {} <= Position <= {}'.format(start, end)
        result = microhapdb.variants.query(query)
        if len(result) > 0:
            yield result


def id_xref(idvalue):
    return next(fetch_by_id(idvalue))
