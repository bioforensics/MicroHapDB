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

import microhapdb
from microhapdb import Population
from textwrap import dedent


def main(args):
    num_markers = len(microhapdb.markers)
    markers = microhapdb.Marker.objectify(microhapdb.markers)
    num_loci = len(set([m.locus for m in markers]))
    print("[microhaplotypes]")
    print(f"  - {num_markers} marker defintions")
    print(f"  - {num_loci} distinct loci")
    num_pops = len(microhapdb.populations)
    num_haplotypes = len(microhapdb.frequencies.groupby(["Marker", "Allele"]))
    num_frequencies = len(microhapdb.frequencies)
    print("[frequencies]")
    print(f"  - {num_haplotypes} haplotypes")
    print(f"  - {num_pops} population groups")
    print(f"  - {num_frequencies} total microhap frequencies")


def subparser(subparsers):
    subparser = subparsers.add_parser(
        "summarize",
        description="Summarize MicroHapDB database contents",
    )
