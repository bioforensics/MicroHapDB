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
    """Retrieve data with a Pandas query.

    Specify the table to query (variant, locus, population, allele) and provide
    a Pandas-style query.

    >>> query_mode('locus', 'Chrom == "chr8"')
                  ID Reference Chrom     Start       End  Source
    77   MHDBL000078    GRCh38  chr8  11738320  11738460  ALFRED
    140  MHDBL000141    GRCh38  chr8   3659270   3659482  ALFRED
    >>> querystr = 'Locus == "MHDBL000177" and Population == "MHDBP000092"'
    >>> query_mode('allele', querystr)
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
    if len(result) == 0:
        return
    print(result.to_string())


def id_mode(idstr):
    """Retrieve data using internal or external IDs, names, or labels.

    >>> id_mode('SI605775E')
                      ID Reference  Chrom  Position Alleles    Source
    9961  MHDBV000009962    GRCh38  chr13  50313423     C,T  dbSNP151
    >>> id_mode('rs690302')
                       ID Reference  Chrom  Position  Alleles    Source
    17205  MHDBV000017206    GRCh38  chr18   8892896  A,C,G,T  dbSNP151
    >>> id_mode('mh19CP-007')
                 ID Reference  Chrom     Start       End  Source
    66  MHDBL000067    GRCh38  chr19  14310740  14310781  ALFRED
    >>> id_mode('SA004109O')
                 ID       Name  Source
    76  MHDBP000077  Colombian  ALFRED
    """
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
    if len(result) == 0:
        return
    print(result.to_string())


def region_mode(region, table=None):
    """Retrieve microhap loci and proximal variants with range queries.

    Use queries of the format "chrX" or "chrX:YYYY-ZZZZ" to retrieve loci or
    variants (or both) from the specified genomic region.

    MicroHapDB includes not only the dbSNP variants that define each
    microhaplotype locus, but also the other variants within its extent and the
    flanking nucleotides.

    >>> region_mode('chr13', table='locus')
                  ID Reference  Chrom      Start        End  Source
    44   MHDBL000045    GRCh38  chr13   24343963   24343994  ALFRED
    50   MHDBL000051    GRCh38  chr13   23191402   23191542  ALFRED
    51   MHDBL000052    GRCh38  chr13   46291795   46291987  ALFRED
    52   MHDBL000053    GRCh38  chr13   53486692   53486837  ALFRED
    68   MHDBL000069    GRCh38  chr13  110154352  110154505  ALFRED
    69   MHDBL000070    GRCh38  chr13   66138600   66138696  ALFRED
    70   MHDBL000071    GRCh38  chr13   94894396   94894513  ALFRED
    116  MHDBL000117    GRCh38  chr13   50313424   50313589  ALFRED
    >>> region_mode('chr18:24557355-24557490')
                 ID Reference  Chrom     Start       End  Source
    34  MHDBL000035    GRCh38  chr18  24557355  24557490  ALFRED
                       ID Reference  Chrom  Position Alleles    Source
    17343  MHDBV000017344    GRCh38  chr18  24557371     C,T  dbSNP151
    17344  MHDBV000017345    GRCh38  chr18  24557373     C,T  dbSNP151
    17345  MHDBV000017346    GRCh38  chr18  24557375     C,T  dbSNP151
    17346  MHDBV000017347    GRCh38  chr18  24557378     C,G  dbSNP151
    17347  MHDBV000017348    GRCh38  chr18  24557389     G,T  dbSNP151
    17348  MHDBV000017349    GRCh38  chr18  24557395     A,G  dbSNP151
    17349  MHDBV000017350    GRCh38  chr18  24557400     C,T  dbSNP151
    17350  MHDBV000017351    GRCh38  chr18  24557402     A,C  dbSNP151
    17351  MHDBV000017352    GRCh38  chr18  24557414     A,T  dbSNP151
    17352  MHDBV000017353    GRCh38  chr18  24557416     A,G  dbSNP151
    17353  MHDBV000017354    GRCh38  chr18  24557431     A,G  dbSNP151
    17354  MHDBV000017355    GRCh38  chr18  24557443     A,G  dbSNP151
    17355  MHDBV000017356    GRCh38  chr18  24557447   C,G,T  dbSNP151
    17356  MHDBV000017357    GRCh38  chr18  24557448     A,G  dbSNP151
    17357  MHDBV000017358    GRCh38  chr18  24557460     C,T  dbSNP151
    17358  MHDBV000017359    GRCh38  chr18  24557463     A,T  dbSNP151
    17359  MHDBV000017360    GRCh38  chr18  24557477     C,T  dbSNP151
    17360  MHDBV000017361    GRCh38  chr18  24557483     A,G  dbSNP151
    17361  MHDBV000017362    GRCh38  chr18  24557486     A,G  dbSNP151
    17362  MHDBV000017363    GRCh38  chr18  24557488     C,T  dbSNP151
    17363  MHDBV000017364    GRCh38  chr18  24557489     G,T  dbSNP151
    """
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
        if len(result) != 0:
            print(result.to_string())
    if table in (None, 'variant'):
        query = 'Chrom == "{}"'.format(chrom)
        if start is not None:
            query += ' and {} < Position < {}'.format(start, end)
        result = microhapdb.variants.query(query)
        if len(result) != 0:
            print(result.to_string())
