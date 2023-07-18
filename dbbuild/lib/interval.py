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

from collections import defaultdict
from intervaltree import IntervalTree


class IntervalIndex:
    def __init__(self):
        self.trees = defaultdict(IntervalTree)
        self.mergeables = dict()

    def add(self, marker):
        self.trees[marker.chrom][marker.start-1:marker.end] = marker

    def check(self):
        for chrom, tree in sorted(self.trees.items()):
            tree.merge_overlaps(data_reducer=lambda locus, marker: locus + [marker], data_initializer=list())
            for interval in tree:
                loci = set([m.locus for m in interval.data])
                if len(loci) == 1:
                    continue
                markers = sorted(interval.data, key=lambda m: m.sortkey)
                for marker in markers[1:]:
                    if marker.name != markers[0].name and marker.name not in self.mergeables:
                        self.mergeables[marker.name] = markers[0].name
                        print(f"Merging '{marker.name}' --> '{markers[0].name}'")
