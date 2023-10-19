#!/usr/bin/env python
# -------------------------------------------------------------------------------------------------
# Copyright (c) 2023, DHS.
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
from collections import defaultdict
import microhapdb
from multiprocessing import Pool
from os import cpu_count
import pandas as pd
import sys


def main(rmsk_path, marker_path, cpus=0, delta=25):
    repeats = parse_ucsc_rmsk_track(rmsk_path)
    marker_table = pd.read_csv(marker_path)
    markers = markers_by_chrom(table=marker_table)
    if cpus < 1:
        cpus = cpu_count()
    data = list(bind_data(repeats, markers, delta=delta))
    with Pool(cpus) as p:
        flagged = p.map(flag_markers, data)
    outdata = list()
    for flagged_markers in flagged:
        outdata.extend(sorted(flagged_markers))
    outframe = pd.DataFrame(outdata, columns=["Marker", "Repeat"])
    outframe = outframe.sort_values(["Marker"])
    return outframe


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
    table = pd.read_csv(path, sep="\t", names=header)
    return table.groupby("genoName")


def markers_by_chrom(table=None):
    if table is None:
        table = microhapdb.markers
    markers = microhapdb.Marker.objectify(table)
    marker_map = defaultdict(list)
    for marker in markers:
        marker_map[marker.chrom].append(marker)
    return marker_map


def bind_data(repeats, markers, delta=0):
    for chrom, repeat_list in repeats:
        if chrom not in markers:
            continue
        marker_list = markers[chrom]
        yield repeat_list, marker_list, delta


def flag_markers(pdata):
    repeats, markers, delta = pdata
    flagged = list()
    for marker in markers:
        start = marker.start - delta
        end = marker.end + delta
        overlap = repeats[(repeats.genoEnd > start) & (repeats.genoStart < end)]
        if len(overlap) > 0:
            marker_repeats = list()
            for i, repeat in overlap.iterrows():
                annot = f"{repeat.repName}|{repeat.repClass}"
                if repeat.repFamily != repeat.repClass and repeat.repFamily != repeat.repName:
                    annot += f"|{repeat.repFamily}"
                marker_repeats.append(annot)
            annots = ";".join(marker_repeats)
            flagged.append((marker.name, annots))
    return flagged


def get_parser():
    parser = ArgumentParser()
    parser.add_argument("rmsk", help="RepeatMasker track from UCSC genome browser (rmsk.txt.gz)")
    parser.add_argument("markers", help="table of markers")
    parser.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        default=sys.stdout,
        help="write output to FILE; by default output is printed to terminal (standard output)",
    )
    parser.add_argument(
        "-c",
        "--cpus",
        type=int,
        default=0,
        metavar="C",
        help="CPUs to use for multiprocessing; by default all available CPUs are used",
    )
    parser.add_argument(
        "-d",
        "--delta",
        type=int,
        default=25,
        metavar="D",
        help="extend D bp beyond the first and last SNPs of the marker when testing for overlap with repeats; by default D=25",
    )
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    table = main(args.rmsk, args.markers, cpus=args.cpus)
    table.to_csv(args.out, index=False)
