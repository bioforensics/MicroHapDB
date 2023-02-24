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

from .tables import markers, populations, frequencies, indels, variantmap, hg38
from .population import Population
from .marker import Marker
from microhapdb import cli
from microhapdb import panel
import pandas as pd
from pkg_resources import resource_filename
from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions


def data_file(path):
    return resource_filename("microhapdb", f"data/{path}")


def set_ae_population(popid="1KGP"):
    global markers
    markers["Ae"] = None
    aes = pd.read_csv(data_file("marker-aes.csv"))
    if popid not in aes.Population.unique():
        raise ValueError(f'no Ae data for population "{popid}"')
    popaes = aes[aes.Population == popid].drop(columns=["Population"])
    markers = markers.drop(columns=["Ae"]).join(popaes.set_index("Marker"), on="Name")


def retrieve_by_id(ident):
    """Retrieve records by name or identifier

    >>> retrieve_by_id("mh17KK-014")
               Name  NumVars  Extent  Chrom    Start      End                Positions              Positions37                          RSIDs    Source     Ae
    609  mh17KK-014        3      37  chr17  4497061  4497097  4497061;4497089;4497097  4400356;4400384;4400392  rs333113;rs8074965;rs11657785  Kidd2018  3.914
    >>> retrieve_by_id("rs8074965")
               Name  NumVars  Extent  Chrom    Start      End                Positions              Positions37                          RSIDs    Source     Ae
    609  mh17KK-014        3      37  chr17  4497061  4497097  4497061;4497089;4497097  4400356;4400384;4400392  rs333113;rs8074965;rs11657785  Kidd2018  3.914
    >>> retrieve_by_id("Chagga")
                   ID    Name        Source
    16  mMHseq-Chagga  Chagga  Gandotra2020
    17      SA000487T  Chagga      Kidd2018
    >>> retrieve_by_id("Asia")
                     ID  Name          Source
    9  MHDBP-936bc36f79  Asia  vanderGaag2018
    >>> retrieve_by_id("Japanese")
                      ID      Name       Source
    54         SA000010B  Japanese     Kidd2018
    55  MHDBP-63967b883e  Japanese  Hiroaki2015
    """

    def id_in_series(ident, series):
        return series.str.contains(ident).any()

    id_in_pop_ids = id_in_series(ident, populations.ID)
    id_in_pop_names = id_in_series(ident, populations.Name)
    id_in_variants = id_in_series(ident, variantmap.Variant)
    id_in_marker_names = id_in_series(ident, markers.Name)
    if id_in_pop_ids or id_in_pop_names:
        return Population.table_from_ids([ident])
    elif id_in_variants or id_in_marker_names:
        return Marker.table_from_ids([ident])
    else:
        raise ValueError(f'identifier "{ident}" not found in MicroHapDB')


set_ae_population("1KGP")
