#!/usr/bin/env python

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
from enum import Enum
import pandas as pd
from util import load_markers

MAX_SIZE = 1.0

# See https://github.com/jordanlab/tagore#input-file-description
def main(markers, marker_aes, markers_passed_filter, markers_failed_filter, final_panel):
    markers = load_markers(markers, marker_aes, objectify=False)

    failed = pd.read_csv(markers_failed_filter, sep="\t")
    failed_markers = markers[markers.Name.isin(failed.Marker)]
    for i, row in failed_markers.iterrows():
        print(row.Chrom, row.Start, row.End, Shape.LINE.value, MAX_SIZE, "#990000", 2, sep="\t")

    passed = pd.read_csv(markers_passed_filter)
    for i, row in passed.iterrows():
        print(row.Chrom, row.Start, row.End, Shape.LINE.value, MAX_SIZE, "#009900", 1, sep="\t")

    final = pd.read_csv(final_panel, sep="\t")
    final_panel = markers[(markers.Name.isin(final.Marker)) | (markers.Chrom == "chrX")]
    max_ae = final_panel.Ae.max()
    for i, row in final_panel.iterrows():
        size = 0.8 * row.Ae / max_ae
        print(row.Chrom, row.Start, row.End, Shape.CIRCLE.value, size, "#0000ff", 1, sep="\t")


class Shape(Enum):
    RECTANGLE = 0
    CIRCLE = 1
    TRIANGLE = 2
    LINE = 3


def get_parser():
    parser = ArgumentParser()
    parser.add_argument("markers", help="path to MicroHapDB marker definitions in CSV format")
    parser.add_argument("aes", help="path to MicroHapDB Ae table in CSV format")
    parser.add_argument("passed", help="path to file of markers that passed filters")
    parser.add_argument("failed", help="path to file of markers that failed filters")
    parser.add_argument("panel", help="path to file containing final panel")
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args.markers, args.aes, args.passed, args.failed, args.panel)
