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


def load_markers(marker_path, ae_path):
    mtable = pd.read_csv(marker_path)
    mtable["Ae"] = 0.0
    if ae_path:
        aes = pd.read_csv(ae_path)
        popaes = aes[aes.Population == "1KGP"].drop(columns=["Population"])
        mtable = mtable.drop(columns=["Ae"]).join(popaes.set_index("Marker"), on="Name")
    markers = sorted(microhapdb.Marker.objectify(mtable), key=lambda m: m.name)
    return markers


def parse_ucsc_rmsk_track(path, filter_key=None, filter_value=None):
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
    return pl.read_csv(path, sep="\t", new_columns=header, has_header=False)
