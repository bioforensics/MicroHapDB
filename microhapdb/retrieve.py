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
                  ID Reference Chrom     Start       End   AvgAe                            Source
    238  MHDBM000239    GRCh38  chr8   1194352   1194372  2.6273  doi:10.1016/j.fsigen.2018.05.008
    239  MHDBM000240    GRCh38  chr8   3659269   3659482  3.9708                            ALFRED
    240  MHDBM000241    GRCh38  chr8  11738319  11738460  2.2865                            ALFRED
    >>> querystr = 'Marker == "MHDBM000213" and Population == "MHDBP000093"'
    >>> for result in fetch_by_query('allele', querystr): print(result)
                Marker   Population Allele  Frequency
    66015  MHDBM000213  MHDBP000093    C,C      0.661
    66016  MHDBM000213  MHDBP000093    C,T      0.339
    66017  MHDBM000213  MHDBP000093    T,C      0.000
    66018  MHDBM000213  MHDBP000093    T,T      0.000
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
    10265  MHDBV000010266    GRCh38  chr13  50313423     C,T  dbSNP151
    >>> for result in fetch_by_id('rs690302'): print(result)
                       ID Reference  Chrom  Position  Alleles    Source
    18914  MHDBV000018915    GRCh38  chr18   8892896  A,C,G,T  dbSNP151
    >>> for result in fetch_by_id('mh06PK-25713'): print(result)
                  ID Reference Chrom     Start       End   AvgAe                            Source
    228  MHDBM000229    GRCh38  chr6  31196949  31197002  3.0005  doi:10.1016/j.fsigen.2018.05.008
    >>> for result in fetch_by_id('mh19CP-007'): print(result)
                  ID Reference  Chrom     Start       End   AvgAe  Source
    130  MHDBM000131    GRCh38  chr19  14310739  14310781  3.0813  ALFRED
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
                 ID Reference  Chrom      Start        End   AvgAe         Source
    65  MHDBM000066    GRCh38  chr13   23191401   23191542  3.6440         ALFRED
    66  MHDBM000067    GRCh38  chr13   23191461   23191542  2.5124  ISFG2019:P597
    67  MHDBM000068    GRCh38  chr13   24343962   24343994  3.0655         ALFRED
    68  MHDBM000069    GRCh38  chr13   46291794   46291987  4.0035         ALFRED
    69  MHDBM000070    GRCh38  chr13   46291948   46291987  2.6806  ISFG2019:P597
    70  MHDBM000071    GRCh38  chr13   50313423   50313589  2.3932         ALFRED
    71  MHDBM000072    GRCh38  chr13   53486691   53486837  6.0444         ALFRED
    72  MHDBM000073    GRCh38  chr13   66138599   66138696  3.4361         ALFRED
    73  MHDBM000074    GRCh38  chr13   94894395   94894513  2.0150         ALFRED
    74  MHDBM000075    GRCh38  chr13  110154351  110154412  3.2411  ISFG2019:P597
    75  MHDBM000076    GRCh38  chr13  110154351  110154505  3.7978         ALFRED
    >>> for result in fetch_by_region('chr18:24557400-24557450'): print(result)
                  ID Reference  Chrom     Start       End   AvgAe         Source
    124  MHDBM000125    GRCh38  chr18  24557354  24557490  2.6580         ALFRED
    125  MHDBM000126    GRCh38  chr18  24557431  24557490  3.1874  ISFG2019:P597
                       ID Reference  Chrom  Position Alleles    Source
    19058  MHDBV000019059    GRCh38  chr18  24557400     C,T  dbSNP151
    19059  MHDBV000019060    GRCh38  chr18  24557402     A,C  dbSNP151
    19060  MHDBV000019061    GRCh38  chr18  24557414     A,T  dbSNP151
    19061  MHDBV000019062    GRCh38  chr18  24557416     A,G  dbSNP151
    19062  MHDBV000019063    GRCh38  chr18  24557431     A,G  dbSNP151
    19063  MHDBV000019064    GRCh38  chr18  24557443     A,G  dbSNP151
    19064  MHDBV000019065    GRCh38  chr18  24557447   C,G,T  dbSNP151
    19065  MHDBV000019066    GRCh38  chr18  24557448     A,G  dbSNP151
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

    >>> standardize_ids(['MHDBM000114'])
    ['MHDBM000114']
    >>> standardize_ids(['SI664623B'])
    ['MHDBM000114']
    >>> standardize_ids(['rs547950691', 'mh02KK-131', 'SA002765U'])
    ['MHDBM000140', 'MHDBP000054', 'MHDBV000027667']
    """
    idmap = microhapdb.idmap
    iddata = idmap[(idmap.XRef.isin(idlist)) | (idmap.mhdbID.isin(idlist))]
    if len(iddata) == 0:
        return list()
    return sorted(iddata.mhdbID.unique())
