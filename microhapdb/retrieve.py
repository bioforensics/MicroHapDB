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
                 Name          PermID Reference  Chrom                                            Offsets      Ae      In     Fst                          Source
    458  mh19USC-19pA  MHDBM-2d713eab    GRCh38  chr19                               561778,561795,561815  2.7453  0.0870  0.0733    10.1016/j.fsigen.2019.102213
    459    mh19KK-056  MHDBM-d6ff8635    GRCh38  chr19                                    4852124,4852324  2.4391  0.0773  0.1760                          ALFRED
    460    mh19CP-007  MHDBM-49dbcc57    GRCh38  chr19                         14310739,14310772,14310780  3.0813  0.0466  0.0776                          ALFRED
    461  mh19USC-19pB  MHDBM-76427274    GRCh38  chr19       16040864,16040894,16040899,16040921,16040929  3.5107  0.1647  0.0731    10.1016/j.fsigen.2019.102213
    462   mh19ZHA-006  MHDBM-5f597339    GRCh38  chr19  20579862,20579863,20579879,20579892,20579916,2...  3.8739  0.1681 -0.0032    10.1016/j.fsigen.2020.102255
    463     mh19NH-23  MHDBM-dd72537b    GRCh38  chr19                         22052723,22052774,22052817  1.9380  0.0414 -0.0124  10.1016/j.legalmed.2015.06.003
    464  mh19KKCS-299  MHDBM-a70896aa    GRCh38  chr19  22546697,22546702,22546748,22546779,22546810,2...     NaN     NaN     NaN    10.1016/j.fsigen.2020.102275
    465    mh19KK-299  MHDBM-8cbeb11c    GRCh38  chr19       22546697,22546748,22546779,22546810,22546850  4.1592  0.2335  0.1102                          ALFRED
    466     mh19AT-47  MHDBM-8f439540    GRCh38  chr19                         22546697,22546748,22546779  2.4025  0.1298  0.1170                   ISFG2019:P597
    467  mh19USC-19qA  MHDBM-f757e745    GRCh38  chr19                33273771,33273785,33273811,33273816  3.3219  0.1127  0.0880    10.1016/j.fsigen.2019.102213
    468    mh19KK-301  MHDBM-2069446a    GRCh38  chr19                50938487,50938502,50938526,50938550  1.9707  0.2673  0.1698                          ALFRED
    469  mh19KKCS-300  MHDBM-bc8b7213    GRCh38  chr19  50947786,50947789,50947790,50947830,50947876,5...     NaN     NaN     NaN    10.1016/j.fsigen.2020.102275
    470    mh19KK-057  MHDBM-eb558c37    GRCh38  chr19                         51654948,51655025,51655062  2.3266  0.0667  0.0428                          ALFRED
    471  mh19USC-19qB  MHDBM-7b40359b    GRCh38  chr19                         53714387,53714389,53714413  4.0756  0.1368  0.0163    10.1016/j.fsigen.2019.102213
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
               Name          PermID Reference  Chrom                  Offsets      Ae      In     Fst  Source
    440  mh17KK-014  MHDBM-83a239de    GRCh38  chr17  4497060,4497088,4497096  2.0215  0.6423  0.3014  ALFRED
    >>> by_id('SI664726F')
               Name          PermID Reference  Chrom                  Offsets      Ae      In     Fst  Source
    440  mh17KK-014  MHDBM-83a239de    GRCh38  chr17  4497060,4497088,4497096  2.0215  0.6423  0.3014  ALFRED
    >>> by_id('MHDBM-ea520d26')
               Name          PermID Reference  Chrom                              Offsets      Ae      In     Fst  Source
    448  mh18KK-285  MHDBM-ea520d26    GRCh38  chr18  24557354,24557431,24557447,24557489  2.7524  0.1721  0.0836  ALFRED
    >>> by_id('PJL')
         ID                           Name Source
    82  PJL  Punjabi from Lahore, Pakistan   1KGP
    >>> by_id('Asia')
                     ID  Name                        Source
    7  MHDBP-936bc36f79  Asia  10.1016/j.fsigen.2018.05.008
    >>> by_id('Japanese')
                      ID      Name                          Source
    45  MHDBP-63967b883e  Japanese  10.1016/j.legalmed.2015.06.003
    46         SA000010B  Japanese                          ALFRED
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
