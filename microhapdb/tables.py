# -------------------------------------------------------------------------------------------------
# Copyright (c) 2022, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from pkg_resources import resource_filename
import pandas as pd


def read_table(table_path):
    table = pd.read_csv(resource_filename("microhapdb", f"data/{table_path}"))
    return table


def compile_variant_map(markers):
    variantmap = markers[["RSIDs", "Name"]].copy()
    variantmap = variantmap.rename(columns={"RSIDs": "Variant", "Name": "Marker"})
    variantmap["Variant"] = variantmap["Variant"].apply(lambda x: x.split(";"))
    return variantmap.explode("Variant")


markers = read_table("marker.csv")
populations = read_table("population.csv")
frequencies = read_table("frequency.csv")
indels = read_table("indels.csv")
variantmap = compile_variant_map(markers)
