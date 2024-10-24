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

from argparse import ArgumentParser
import microhapdb
import pandas as pd
import polars as pl
import sys
from tqdm import tqdm


def main(markers, repeats):
    data = list()
    for mh in tqdm(markers):
        for offset in mh.offsets:
            start = offset - 4
            end = offset + 4
            overlapping = repeats.filter(
                (pl.col("genoName") == mh.chrom)
                & (pl.col("genoStart") < end)
                & (start < pl.col("genoEnd"))
            )
            if len(overlapping) > 0:
                repeat_data = "|".join(overlapping.get_column("repName"))
                entry = (mh.name, len(mh), mh.data.Ae, offset, repeat_data)
                data.append(entry)
    output = pd.DataFrame(data, columns=["Marker", "Extent", "Ae", "Offset", "Repeats"])
    output.to_csv(sys.stdout, sep="\t", index=False)


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
    return pl.read_csv(path, sep="\t", new_columns=header, has_header=False)


def get_parser():
    parser = ArgumentParser()
    parser.add_argument("markers")
    parser.add_argument("rmsk_path")
    parser.add_argument("--aes")
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    marker_table = pd.read_csv(args.markers)
    marker_table["Ae"] = 0.0
    if args.aes:
        aes = pd.read_csv(args.aes)
        popaes = aes[aes.Population == "1KGP"].drop(columns=["Population"])
        marker_table = marker_table.drop(columns=["Ae"]).join(
            popaes.set_index("Marker"), on="Name"
        )
    markers = list(microhapdb.Marker.objectify(marker_table))
    repeats = parse_ucsc_rmsk_track(args.rmsk_path)
    repeats = repeats.filter(pl.col("repClass") == "Low_complexity")
    main(markers, repeats)
