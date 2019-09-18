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

    >>> for result in fetch_by_query('marker', 'Chrom == "chr8"'): print(result)
                  ID Reference Chrom     Start       End   AvgAe  Source
    202  MHDBM000203    GRCh38  chr8   1194352   1194372  2.6273    LOVD
    203  MHDBM000204    GRCh38  chr8   3659269   3659482  3.9708  ALFRED
    204  MHDBM000205    GRCh38  chr8  11738319  11738460  2.2865  ALFRED
    >>> querystr = 'Marker == "MHDBM000177" and Population == "MHDBP000092"'
    >>> for result in fetch_by_query('allele', querystr): print(result)
                Marker   Population Allele  Frequency
    66015  MHDBM000177  MHDBP000092    C,C      0.661
    66016  MHDBM000177  MHDBP000092    C,T      0.339
    66017  MHDBM000177  MHDBP000092    T,C      0.000
    66018  MHDBM000177  MHDBP000092    T,T      0.000
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
    10491  MHDBV000010492    GRCh38  chr13  50313423     C,T  dbSNP151
    >>> for result in fetch_by_id('rs690302'): print(result)
                       ID Reference  Chrom  Position  Alleles    Source
    19313  MHDBV000019314    GRCh38  chr18   8892896  A,C,G,T  dbSNP151
    >>> for result in fetch_by_id('mh06PK-25713'): print(result)
                  ID Reference Chrom     Start       End   AvgAe Source
    192  MHDBM000193    GRCh38  chr6  31196949  31197002  3.0005   LOVD
    >>> for result in fetch_by_id('mh19CP-007'): print(result)
                  ID Reference  Chrom     Start       End   AvgAe  Source
    109  MHDBM000110    GRCh38  chr19  14310739  14310781  3.0813  ALFRED
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
    markers.

    >>> for result in fetch_by_region('chr13', table='marker'): print(result)
                 ID Reference  Chrom      Start        End   AvgAe  Source
    55  MHDBM000056    GRCh38  chr13   23191401   23191542  3.6440  ALFRED
    56  MHDBM000057    GRCh38  chr13   24343962   24343994  3.0655  ALFRED
    57  MHDBM000058    GRCh38  chr13   46291794   46291987  4.0035  ALFRED
    58  MHDBM000059    GRCh38  chr13   50313423   50313589  2.3932  ALFRED
    59  MHDBM000060    GRCh38  chr13   53486691   53486837  6.0444  ALFRED
    60  MHDBM000061    GRCh38  chr13   66138599   66138696  3.4361  ALFRED
    61  MHDBM000062    GRCh38  chr13   94894395   94894513  2.0150  ALFRED
    62  MHDBM000063    GRCh38  chr13  110154351  110154505  3.7978  ALFRED
    >>> for result in fetch_by_region('chr18:24557400-24557450'): print(result)
                  ID Reference  Chrom     Start       End  AvgAe  Source
    105  MHDBM000106    GRCh38  chr18  24557354  24557490  2.658  ALFRED
                       ID Reference  Chrom  Position Alleles    Source
    19457  MHDBV000019458    GRCh38  chr18  24557400     C,T  dbSNP151
    19458  MHDBV000019459    GRCh38  chr18  24557402     A,C  dbSNP151
    19459  MHDBV000019460    GRCh38  chr18  24557414     A,T  dbSNP151
    19460  MHDBV000019461    GRCh38  chr18  24557416     A,G  dbSNP151
    19461  MHDBV000019462    GRCh38  chr18  24557431     A,G  dbSNP151
    19462  MHDBV000019463    GRCh38  chr18  24557443     A,G  dbSNP151
    19463  MHDBV000019464    GRCh38  chr18  24557447   C,G,T  dbSNP151
    19464  MHDBV000019465    GRCh38  chr18  24557448     A,G  dbSNP151
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
    if table not in (None, 'marker', 'variant'):
        msg = 'region query not supported for table "{}"'.format(table)
        raise ValueError(msg)
    if table in (None, 'marker'):
        query = 'Chrom == "{}"'.format(chrom)
        if start is not None:
            query += ' and (Start < {}'.format(end)
            query += ' and End > {})'.format(start)
        result = microhapdb.markers.query(query)
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


def allele_positions(marker):
    markers = microhapdb.id_xref(marker)
    assert len(markers) == 1
    markerid = markers.iloc[0].ID
    variants = microhapdb.variants[
        microhapdb.variants.ID.isin(
            microhapdb.variantmap[
                microhapdb.variantmap.MarkerID == markerid
            ].VariantID
        )
    ]
    return list(variants.drop_duplicates('Position').Position)


def standardize_ids(idlist):
    """Convert a list of IDs or DB cross-references to internal MicroHapDB IDs.

    >>> standardize_ids(['MHDBM000097'])
    ['MHDBM000097']
    >>> standardize_ids(['SI664623B'])
    ['MHDBM000097']
    >>> standardize_ids(['rs547950691', 'mh02KK-131', 'SA002765U'])
    ['MHDBM000118', 'MHDBP000054', 'MHDBV000028356']
    """
    idmap = microhapdb.idmap
    iddata = idmap[(idmap.XRef.isin(idlist)) | (idmap.mhdbID.isin(idlist))]
    if len(iddata) == 0:
        return list()
    return sorted(iddata.mhdbID.unique())
