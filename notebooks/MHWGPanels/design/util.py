# -------------------------------------------------------------------------------------------------
# Copyright (c) 2024, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

import microhapdb
import pandas as pd
import polars as pl


def load_markers(marker_path, ae_path, objectify=True):
    mtable = pd.read_csv(marker_path)
    mtable["Ae"] = 0.0
    if ae_path:
        aes = pd.read_csv(ae_path)
        popaes = aes[aes.Population == "1KGP"].drop(columns=["Population"])
        mtable = mtable.drop(columns=["Ae"]).join(popaes.set_index("Marker"), on="Name")
    if objectify:
        return sorted(microhapdb.Marker.objectify(mtable), key=lambda m: m.name)
    else:
        return mtable


def parse_ucsc_rmsk_track(path):
    header = [
        "bin",
        "swScore",
        "milliDiv",
        "milliDel",
        "milliIns",
        "genoName",
        "genoStart",
        "genoEnd",
        "genoLeft",
        "strand",
        "repName",
        "repClass",
        "repFamily",
        "repStart",
        "repEnd",
        "repLeft",
        "id",
    ]
    return pl.read_csv(path, separator="\t", new_columns=header, has_header=False)


def marker_distance(m1, m2):
    return interval_distance((m1.start, m1.end), (m2.start, m2.end))


def interval_distance(itvl1, itvl2):
    x, y = sorted((itvl1, itvl2))
    distance = y[0] - x[1] if x[0] <= x[1] < y[0] else 0
    return distance
