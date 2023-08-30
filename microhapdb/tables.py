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
from pathlib import Path
from pyfaidx import Fasta as FastaIdx


def read_table(table_path):
    table = pd.read_csv(resource_filename("microhapdb", f"data/{table_path}"))
    return table


def compile_variant_map(markers):
    variantmap = markers[["RSIDs", "Name"]].copy()
    variantmap = variantmap.rename(columns={"RSIDs": "Variant", "Name": "Marker"})
    variantmap["Variant"] = variantmap["Variant"].apply(
        lambda x: [] if pd.isna(x) else x.split(";")
    )
    return variantmap.explode("Variant")


markers = read_table("marker.csv")
markers["Ae"] = None
merged = read_table("merged.csv")
populations = read_table("population.csv")
frequencies = read_table("frequency.csv.gz")
indels = read_table("indels.csv")
variantmap = compile_variant_map(markers)
hg38file = resource_filename("microhapdb", "data/hg38.fasta")
if Path(hg38file).is_file():
    hg38 = FastaIdx(resource_filename("microhapdb", "data/hg38.fasta"))
else:  # pragma: no cover
    hg38 = None
