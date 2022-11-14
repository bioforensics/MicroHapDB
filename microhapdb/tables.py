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
    table = pd.read_csv(resource_filename("microhapdb", f"data/{table_path}"), sep="\t")
    return table


markers = read_table("marker.tsv")
populations = read_table("population.tsv")
frequencies = read_table("frequency.tsv")
variantmap = read_table("variantmap.tsv")
idmap = read_table("idmap.tsv")
sequences = read_table("sequences.tsv")
indels = read_table("indels.tsv")
