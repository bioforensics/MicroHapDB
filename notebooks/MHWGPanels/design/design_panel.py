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
from linkage_graph import LinkageGraph, LinkageGraphThresholds, logger
from util import load_markers, marker_distance


def main(markers, thresholds, cutlist=None):
    logger.info("Populating linkage graph")
    graph = LinkageGraph.populate(markers, thresholds)
    panel = graph.design_panel()
    print("Chromosome", "Marker", "Extent", "Ae", sep="\t")
    for mh in panel:
        print(mh.chrom, mh.name, len(mh), mh.data.Ae, sep="\t")
    if cutlist:
        with open(cutlist, "w") as fh:
            generate_cut_list(fh, panel)


def generate_cut_list(outstream, panel, min_ae=5.0):
    print(
        "Chromosome", "Marker", "Extent", "Ae", "ClosestNeighbor", "NeighborAe", "NeighborDistance",
        sep="\t",
        file=outstream
    )
    panel_loci = set(mh.locus for mh in panel)
    for marker in markers:
        if marker.data.Ae < min_ae or marker.locus in panel_loci:
            continue
        neighbor, distance = find_closest_panel_marker(marker, panel)
        print(
            marker.chrom,
            marker.name,
            len(marker),
            marker.data.Ae,
            neighbor.name,
            neighbor.data.Ae,
            distance,
            sep="\t",
            file=outstream,
        )


def find_closest_panel_marker(marker, panel):
    closest_neighbor = None
    closest_distance = float("Inf")
    for neighbor in panel:
        if neighbor.chrom != marker.chrom:
            continue
        distance = marker_distance(marker, neighbor)
        if distance < closest_distance:
            closest_neighbor = neighbor
            closest_distance = distance
    return closest_neighbor, closest_distance


def get_parser():
    parser = ArgumentParser()
    parser.add_argument("markers", help="path to MicroHapDB marker definitions in CSV format")
    parser.add_argument("aes", help="path to MicroHapDB Ae table in CSV format")
    parser.add_argument(
        "--max-per-chrom",
        type=int,
        nargs=2,
        default=(8, 8),
        metavar="M",
        help="select the M highest ranked short and long markers (respectively) per chromosome by Ae; by default M=(8, 8)",
    )
    parser.add_argument(
        "--batches",
        type=int,
        default=8,
        metavar="B",
        help="split the linkage graph into B batches; by default B=8",
    )
    parser.add_argument(
        "--distance",
        type=float,
        default=10e6,
        metavar="D",
        help="two markers must be separated by more than D bp to be considered independently inherited (as a heuristic); by default D=10000000",
    )
    parser.add_argument(
        "--cut-list",
        metavar="FILE",
        help="write high Ae markers that didn't make the cut to FILE",
    )
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    markers = load_markers(args.markers, args.aes)
    thresholds = LinkageGraphThresholds(
        max_per_chrom_short=args.max_per_chrom[0],
        max_per_chrom_long=args.max_per_chrom[1],
        ld_distance=args.distance,
    )
    main(markers, thresholds, cutlist=args.cut_list)
