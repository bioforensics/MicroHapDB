# -------------------------------------------------------------------------------------------------
# Copyright (c) 2018, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from .tables import markers, populations, frequencies, variantmap, idmap, sequences, indels
from .population import Population
from .marker import Marker
from microhapdb import cli
from microhapdb import panel
import pandas as pd
from pkg_resources import resource_filename
from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions


def set_ae_population(popid=None):
    global markers
    columns = ["Name", "PermID", "Reference", "Chrom", "Offsets", "Ae", "In", "Fst", "Source"]
    if popid is None:
        filename = resource_filename("microhapdb", "data/marker.tsv")
        defaults = pd.read_csv(filename, sep="\t")
        defaults = defaults[["Name", "Ae"]]
        markers = markers.drop(columns=["Ae"]).join(defaults.set_index("Name"), on="Name")[columns]
    else:
        filename = resource_filename("microhapdb", "data/marker-aes.tsv")
        aes = pd.read_csv(filename, sep="\t")
        if popid not in aes.Population.unique():
            raise ValueError(f'no Ae data for population "{popid}"')
        popaes = aes[aes.Population == popid].drop(columns=["Population"])
        markers = markers.drop(columns=["Ae"]).join(popaes.set_index("Marker"), on="Name")[columns]


def set_reference(refr):
    global markers
    assert refr in (37, 38)
    columns = ["Name", "PermID", "Reference", "Chrom", "Offsets", "Ae", "In", "Fst", "Source"]
    if refr == 38:
        filename = resource_filename("microhapdb", "data/marker.tsv")
        defaults = pd.read_csv(filename, sep="\t")[["Name", "Reference", "Offsets"]]
        markers = markers.drop(columns=["Reference", "Offsets"]).join(
            defaults.set_index("Name"), on="Name"
        )[columns]
    else:
        filename = resource_filename("microhapdb", "data/marker-offsets-GRCh37.tsv")
        o37 = pd.read_csv(filename, sep="\t")
        markers = markers.drop(columns=["Reference", "Offsets"]).join(
            o37.set_index("Marker"), on="Name"
        )[columns]


def retrieve_by_id(ident):
    """Retrieve records by name or identifier

    >>> retrieve_by_id("mh17KK-014")
               Name          PermID Reference  Chrom                  Offsets      Ae      In     Fst  Source
    510  mh17KK-014  MHDBM-83a239de    GRCh38  chr17  4497060,4497088,4497096  2.0215  0.6423  0.3014  ALFRED
    >>> retrieve_by_id("SI664726F")
               Name          PermID Reference  Chrom                  Offsets      Ae      In     Fst  Source
    510  mh17KK-014  MHDBM-83a239de    GRCh38  chr17  4497060,4497088,4497096  2.0215  0.6423  0.3014  ALFRED
    >>> retrieve_by_id("MHDBM-ea520d26")
               Name          PermID Reference  Chrom                              Offsets      Ae      In     Fst  Source
    539  mh18KK-285  MHDBM-ea520d26    GRCh38  chr18  24557354,24557431,24557447,24557489  2.7524  0.1721  0.0836  ALFRED
    >>> retrieve_by_id("PJL")
         ID                           Name Source
    82  PJL  Punjabi from Lahore, Pakistan   1KGP
    >>> retrieve_by_id("Asia")
                     ID  Name                        Source
    7  MHDBP-936bc36f79  Asia  10.1016/j.fsigen.2018.05.008
    >>> retrieve_by_id("Japanese")
                      ID      Name                          Source
    45  MHDBP-63967b883e  Japanese  10.1016/j.legalmed.2015.06.003
    46         SA000010B  Japanese                          ALFRED
    """

    def id_in_series(ident, series):
        return series.str.contains(ident).any()

    if id_in_series(ident, idmap.Xref):
        result = idmap[idmap.Xref == ident]
        assert len(result) == 1
        ident = result.ID.iloc[0]
    id_in_pop_ids = id_in_series(ident, populations.ID)
    id_in_pop_names = id_in_series(ident, populations.Name)
    id_in_variants = id_in_series(ident, variantmap.Variant)
    id_in_marker_names = id_in_series(ident, markers.Name)
    id_in_marker_permids = id_in_series(ident, markers.PermID)
    if id_in_pop_ids or id_in_pop_names:
        return Population.table_from_ids([ident])
    elif id_in_variants or id_in_marker_names or id_in_marker_permids:
        return Marker.table_from_ids([ident])
    else:
        raise ValueError(f'identifier "{ident}" not found in MicroHapDB')
