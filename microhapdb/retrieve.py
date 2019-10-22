# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import microhapdb


# ----- [ Supporting internals ] ----- #

def id_in_series(ident, series):
    return series.str.contains(ident).any()


def standardize_marker_ids(idents):
    ids = set()
    for ident in idents:
        if id_in_series(ident, microhapdb.variantmap.Variant):
            markernames = microhapdb.variantmap[microhapdb.variantmap.Variant == ident].Marker
            ids.update(markernames)
        elif id_in_series(ident, microhapdb.markers.PermID):
            result = microhapdb.markers[microhapdb.markers.PermID == ident]
            ids.update(result.Name)
        elif id_in_series(ident, microhapdb.markers.Name):
            result = microhapdb.markers[microhapdb.markers.Name == ident]
            ids.update(result.Name)
        elif id_in_series(ident, microhapdb.idmap.Xref):
            markername = microhapdb.idmap[microhapdb.idmap.Xref == ident].ID
            assert len(markername) == 1
            markername = markername.iloc[0]
            result = microhapdb.markers[microhapdb.markers.Name == markername]
            ids.update(result.Name)
    return microhapdb.markers[microhapdb.markers.Name.isin(ids)].Name


def standardize_population_ids(idents):
    result = microhapdb.populations[
        (microhapdb.populations.ID.isin(idents)) | (microhapdb.populations.Name.isin(idents))
    ]
    return result.ID


# ----- [ Supported API ] ----- #

def by_region(regionstr):
    '''Retrieve data from the specied genomic region

    >>> by_region('chr19')
               Name          PermID Reference  Chrom                                       Offsets   AvgAe                          Source
    256  mh19KK-056  MHDBM-d6ff8635    GRCh38  chr19                               4852124,4852324  2.4143                          ALFRED
    257  mh19CP-007  MHDBM-49dbcc57    GRCh38  chr19                    14310739,14310772,14310780  3.0813                          ALFRED
    258   mh19NH-23  MHDBM-dd72537b    GRCh38  chr19                    22052723,22052774,22052817     NaN  10.1016/j.legalmed.2015.06.003
    259  mh19KK-299  MHDBM-8cbeb11c    GRCh38  chr19  22546697,22546748,22546779,22546810,22546850  3.8989                          ALFRED
    260   mh19AT-47  MHDBM-8f439540    GRCh38  chr19                    22546697,22546748,22546779  1.4537                   ISFG2019:P597
    261  mh19KK-301  MHDBM-2069446a    GRCh38  chr19           50938487,50938502,50938526,50938550  1.8143                          ALFRED
    262  mh19KK-057  MHDBM-eb558c37    GRCh38  chr19                    51654948,51655025,51655062  2.1923                          ALFRED
    >>> by_region('chr18:1-25000000')
                 Name          PermID Reference  Chrom                                          Offsets   AvgAe                        Source
    249  mh18PK-87558  MHDBM-1e5374f1    GRCh38  chr18  1960542,1960557,1960561,1960566,1960582,1960588  2.3650  10.1016/j.fsigen.2018.05.008
    250    mh18CP-005  MHDBM-a85754d3    GRCh38  chr18                  8892864,8892893,8892896,8892907  3.6635                        ALFRED
    251    mh18KK-285  MHDBM-ea520d26    GRCh38  chr18              24557354,24557431,24557447,24557489  2.6580                        ALFRED
    252     mh18AT-38  MHDBM-db09ec41    GRCh38  chr18                       24557431,24557447,24557489  3.1870                 ISFG2019:P597
    '''
    markers = microhapdb.markers.copy()
    markers['Start'] = markers.Offsets.apply(lambda o: min(map(int, o.split(','))))
    markers['End'] = markers.Offsets.apply(lambda o: max(map(int, o.split(','))) + 1)
    chrom, start, end = microhapdb.util.parse_regionstr(regionstr)
    query = 'Chrom == "{}"'.format(chrom)
    if start is not None:
        query += ' and (Start < {}'.format(end)
        query += ' and End > {})'.format(start)
    result = markers.query(query)
    return result.drop(columns=['Start', 'End'])


def by_id(ident):
    '''Retrieve records by name or identifier

    >>> by_id('mh17KK-014')
               Name          PermID Reference  Chrom                  Offsets   AvgAe  Source
    234  mh17KK-014  MHDBM-83a239de    GRCh38  chr17  4497060,4497088,4497096  1.9183  ALFRED
    >>> by_id('SI664726F')
               Name          PermID Reference  Chrom                  Offsets   AvgAe  Source
    234  mh17KK-014  MHDBM-83a239de    GRCh38  chr17  4497060,4497088,4497096  1.9183  ALFRED
    >>> by_id('MHDBM-ea520d26')
               Name          PermID Reference  Chrom                              Offsets  AvgAe  Source
    251  mh18KK-285  MHDBM-ea520d26    GRCh38  chr18  24557354,24557431,24557447,24557489  2.658  ALFRED
    >>> by_id('Asia')
                     ID  Name                        Source
    6  MHDBP-936bc36f79  Asia  10.1016/j.fsigen.2018.05.008
    >>> by_id('Japanese')
                      ID      Name                          Source
    42         SA004060K  Japanese                          ALFRED
    43         SA000010B  Japanese                          ALFRED
    44  MHDBP-63967b883e  Japanese  10.1016/j.legalmed.2015.06.003
    >>> by_id('SA004240K')
               ID     Name  Source
    79  SA004240K  Punjabi  ALFRED
    '''
    if id_in_series(ident, microhapdb.populations.ID) or id_in_series(ident, microhapdb.populations.Name):
        idlist = standardize_population_ids([ident])
        return microhapdb.populations[microhapdb.populations.ID.isin(idlist)]
    series = [
        microhapdb.variantmap.Variant, microhapdb.idmap.Xref, microhapdb.markers.Name,
        microhapdb.markers.PermID
    ]
    for s in series:
        if id_in_series(ident, s):
            idlist = standardize_marker_ids([ident])
            return microhapdb.markers[microhapdb.markers.Name.isin(idlist)]
    else:
        raise ValueError('identifier "{}" not found in MicroHapDB'.format(ident))
