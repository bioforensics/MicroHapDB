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
                  ID Reference Chrom     Start       End  Source
    77   MHDBL000078    GRCh38  chr8  11738320  11738460  ALFRED
    140  MHDBL000141    GRCh38  chr8   3659270   3659482  ALFRED
    >>> querystr = 'Locus == "MHDBL000177" and Population == "MHDBP000092"'
    >>> for result in fetch_by_query('allele', querystr): print(result)
                 Locus   Population   Allele  Frequency
    71439  MHDBL000177  MHDBP000092  T,C,A,A      0.288
    71440  MHDBL000177  MHDBP000092  T,C,G,A      0.025
    71441  MHDBL000177  MHDBP000092  T,C,G,G      0.071
    71442  MHDBL000177  MHDBP000092  T,G,A,A      0.268
    71443  MHDBL000177  MHDBP000092  T,G,A,G      0.005
    71444  MHDBL000177  MHDBP000092  T,G,G,A      0.025
    71445  MHDBL000177  MHDBP000092  T,G,G,G      0.318
    71446  MHDBL000177  MHDBP000092  G,C,A,A      0.000
    71447  MHDBL000177  MHDBP000092  G,C,A,G      0.000
    71448  MHDBL000177  MHDBP000092  G,C,G,A      0.000
    71449  MHDBL000177  MHDBP000092  G,C,G,G      0.000
    71450  MHDBL000177  MHDBP000092  G,G,A,A      0.000
    71451  MHDBL000177  MHDBP000092  G,G,G,A      0.000
    71452  MHDBL000177  MHDBP000092  G,G,G,G      0.000
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
    9961  MHDBV000009962    GRCh38  chr13  50313423     C,T  dbSNP151
    >>> for result in fetch_by_id('rs690302'): print(result)
                       ID Reference  Chrom  Position  Alleles    Source
    17205  MHDBV000017206    GRCh38  chr18   8892896  A,C,G,T  dbSNP151
    >>> for result in fetch_by_id('mh19CP-007'): print(result)
                 ID Reference  Chrom     Start       End  Source
    66  MHDBL000067    GRCh38  chr19  14310740  14310781  ALFRED
    >>> for result in fetch_by_id('SA004109O'): print(result)
                 ID       Name  Source
    76  MHDBP000077  Colombian  ALFRED
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
                  ID Reference  Chrom      Start        End  Source
    44   MHDBL000045    GRCh38  chr13   24343963   24343994  ALFRED
    50   MHDBL000051    GRCh38  chr13   23191402   23191542  ALFRED
    51   MHDBL000052    GRCh38  chr13   46291795   46291987  ALFRED
    52   MHDBL000053    GRCh38  chr13   53486692   53486837  ALFRED
    68   MHDBL000069    GRCh38  chr13  110154352  110154505  ALFRED
    69   MHDBL000070    GRCh38  chr13   66138600   66138696  ALFRED
    70   MHDBL000071    GRCh38  chr13   94894396   94894513  ALFRED
    116  MHDBL000117    GRCh38  chr13   50313424   50313589  ALFRED
    >>> for result in fetch_by_region('chr18:24557400-24557450'): print(result)
                 ID Reference  Chrom     Start       End  Source
    34  MHDBL000035    GRCh38  chr18  24557355  24557490  ALFRED
                       ID Reference  Chrom  Position Alleles    Source
    17349  MHDBV000017350    GRCh38  chr18  24557400     C,T  dbSNP151
    17350  MHDBV000017351    GRCh38  chr18  24557402     A,C  dbSNP151
    17351  MHDBV000017352    GRCh38  chr18  24557414     A,T  dbSNP151
    17352  MHDBV000017353    GRCh38  chr18  24557416     A,G  dbSNP151
    17353  MHDBV000017354    GRCh38  chr18  24557431     A,G  dbSNP151
    17354  MHDBV000017355    GRCh38  chr18  24557443     A,G  dbSNP151
    17355  MHDBV000017356    GRCh38  chr18  24557447   C,G,T  dbSNP151
    17356  MHDBV000017357    GRCh38  chr18  24557448     A,G  dbSNP151
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
