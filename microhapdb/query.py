#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import microhapdb


def query(table, querystr):
    '''Retrieve data using a Pandas-style query.

    >>> query('marker', 'Name.str.contains("PK") and Chrom == "chr6"')
                Name          PermID Reference Chrom                                            Offsets   AvgAe                        Source
    85  mh06PK-24844  MHDBM-aa39cbba    GRCh38  chr6  13861392,13861399,13861414,13861421,13861430,1...  2.4481  10.1016/j.fsigen.2018.05.008
    89  mh06PK-25713  MHDBM-7d00efdc    GRCh38  chr6  31196949,31196961,31196972,31196985,31196992,3...  3.0005  10.1016/j.fsigen.2018.05.008
    >>> query('population', 'Name.str.contains("Afr")')
              ID               Name                        Source
    1     Africa             Africa  10.1016/j.fsigen.2018.05.008
    2  SA004047P  African Americans                        ALFRED
    3  SA000101C  African Americans                        ALFRED
    4  SA004242M    Afro-Caribbeans                        ALFRED
    '''
    if not table:
        raise ValueError('must specify table to invoke a query')
    if table not in microhapdb.tables:
        raise ValueError('unsupported table "{}"'.format(table))
    table = microhapdb.tables[table]
    return table.query(querystr)


def parse_regionstr(regionstr):
    chrom, start, end = None, None, None
    if ':' in regionstr:
        chrom, rng = regionstr.split(':')
        if rng.count('-') != 1:
            raise ValueError('cannot parse region "{}"'.format(regionstr))
        startstr, endstr = rng.split('-')
        start, end = int(startstr), int(endstr)
    else:
        chrom = regionstr
    return chrom, start, end


def region(regionstr):
    '''Retrieve data from the specied genomic region

    >>> region('chr19')
               Name          PermID Reference  Chrom                                       Offsets   AvgAe         Source
    224  mh19KK-056  MHDBM-d6ff8635    GRCh38  chr19                               4852124,4852324  2.4143         ALFRED
    225  mh19CP-007  MHDBM-49dbcc57    GRCh38  chr19                    14310739,14310772,14310780  3.0813         ALFRED
    226  mh19KK-299  MHDBM-8cbeb11c    GRCh38  chr19  22546697,22546748,22546779,22546810,22546850  3.8989         ALFRED
    227   mh19AT-47  MHDBM-8f439540    GRCh38  chr19                    22546697,22546748,22546779  1.4537  ISFG2019:P597
    228  mh19KK-301  MHDBM-2069446a    GRCh38  chr19           50938487,50938502,50938526,50938550  1.8143         ALFRED
    229  mh19KK-057  MHDBM-eb558c37    GRCh38  chr19                    51654948,51655025,51655062  2.1923         ALFRED
    >>> region('chr18:1-25000000')
                 Name          PermID Reference  Chrom                                          Offsets   AvgAe                        Source
    217  mh18PK-87558  MHDBM-1e5374f1    GRCh38  chr18  1960542,1960557,1960561,1960566,1960582,1960588  2.3650  10.1016/j.fsigen.2018.05.008
    218    mh18CP-005  MHDBM-a85754d3    GRCh38  chr18                  8892864,8892893,8892896,8892907  3.6614                        ALFRED
    219    mh18KK-285  MHDBM-ea520d26    GRCh38  chr18              24557354,24557431,24557447,24557489  2.6580                        ALFRED
    220     mh18AT-38  MHDBM-db09ec41    GRCh38  chr18                       24557431,24557447,24557489  3.1870                 ISFG2019:P597
    '''
    markers = microhapdb.markers.copy()
    markers['Start'] = markers.Offsets.apply(lambda o: min(map(int, o.split(','))))
    markers['End'] = markers.Offsets.apply(lambda o: max(map(int, o.split(','))) + 1)
    chrom, start, end = parse_regionstr(regionstr)
    query = 'Chrom == "{}"'.format(chrom)
    if start is not None:
        query += ' and (Start < {}'.format(end)
        query += ' and End > {})'.format(start)
    result = markers.query(query)
    return result.drop(columns=['Start', 'End'])


def id(ident):
    if microhapdb.populations.ID.str.contains(ident).any():
        return microhapdb.populations[microhapdb.populations.ID == ident]
    if microhapdb.variantmap.Variant.str.contains(ident).any():
        markernames = microhapdb.variantmap[microhapdb.variantmap.Variant == ident].Marker
        return microhapdb.markers[microhapdb.markers.Name.isin(markernames)]
    if microhapdb.idmap.Xref.str.contains(ident).any():
        markername = microhapdb.idmap[microhapdb.idmap.Xref == ident].ID
        assert len(markername) == 1
        markername = list(markername)[0]
        return microhapdb.markers[microhapdb.markers.Name == markername]
    if microhapdb.markers.Name.str.contains(ident).any():
        return microhapdb.markers[microhapdb.markers.Name == ident]
    else:
        raise ValueError('identifier "{}" not found in MicroHapDB'.format(ident))
