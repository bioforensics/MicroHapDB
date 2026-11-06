#!/usr/bin/env python

# -------------------------------------------------------------------------------------------------
# Copyright (c) 2026, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

import sys
from argparse import ArgumentParser
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import combinations
from math import comb
from pathlib import Path

from matplotlib import pyplot as plt
from tqdm import tqdm
from util import load_markers


def main(microhaps, outfile):
    distances = DistanceCalculator(microhaps).calculate_distances()
    distances_capped = [d / 1e6 if d < 5e6 else 5 for d in distances]
    fig = plt.figure(figsize=(6, 4), dpi=300)
    plt.hist(distances_capped, bins=40)
    plt.xlabel("Distance (Mb) to closest candidate MH")
    plt.ylabel("Frequency")
    plt.gca().set_axisbelow(True)
    plt.grid(True, linestyle="--", color="#eeeeee")
    fig.savefig(outfile, bbox_inches="tight")


def get_parser():
    parser = ArgumentParser()
    parser.add_argument("microhaps", help="path to final microhap panel in CSV format")
    parser.add_argument("outfile", help="output filename")
    return parser


@dataclass
class DistanceCalculator:
    infile: Path
    markers: defaultdict = field(init=False)

    def __post_init__(self):
        self.infile = Path(self.infile)
        self.markers = defaultdict(lambda: defaultdict(set))
        for marker in load_markers(self.infile, None, objectify=True):
            self.markers[marker.chrom][marker.locus].add(marker)

    def calculate_distances(self):
        distances = {}
        for locus1, locus2 in self.locus_pairs:
            distance = self.pairwise_locus_distance(locus1, locus2)
            for locus in (locus1, locus2):
                locusname = next(iter(locus)).locus
                prev_lowest = distances.get(locusname, float("Inf"))
                distances[locusname] = min(distance, prev_lowest)
        big = [(locus, dist) for locus, dist in distances.items() if dist > 1e7]
        if big:
            print("-----[MH ≥ 10Mb away from any other MH]-----", *big, sep="\n", file=sys.stderr)
        return sorted(distances.values())

    @property
    def locus_pairs(self):
        total_pairs = sum(comb(len(loci), 2) for loci in self.markers.values())
        locus_stream = tqdm(self.chromosome_pairs, total=total_pairs, desc="MH locus pairs")
        yield from locus_stream

    @property
    def chromosome_pairs(self):
        for loci in self.markers.values():
            yield from combinations(loci.values(), 2)

    @staticmethod
    def pairwise_locus_distance(locus1, locus2):
        center1 = DistanceCalculator.locus_midpoint(locus1)
        center2 = DistanceCalculator.locus_midpoint(locus2)
        return abs(center1 - center2)

    @staticmethod
    def locus_midpoint(locus):
        min_coord = min(m.start for m in locus)
        max_coord = max(m.end for m in locus)
        return (min_coord + max_coord) / 2



if __name__ == "__main__":
    args = get_parser().parse_args()
    main(args.microhaps, args.outfile)
