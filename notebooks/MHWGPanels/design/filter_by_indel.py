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
from indel_index import IndelIndex
import sys
from util import load_markers


def main(markers, dbsnp_path, delta=25, min_freq=0.005, distance=4):
    index = IndelIndex()
    index.populate(markers, dbsnp_path, delta=delta, min_freq=min_freq)
    message = f"Found {len(index)} indels with an alt allele at frequency ≥ {min_freq}"
    print(message, file=sys.stderr, flush=True)
    print("Marker", "Extent", "Ae", "MaxAltFreq", "ADSs", "InDels", sep="\t")
    for mh in markers:
        process_marker(mh, index, distance=distance)


def process_marker(mh, index, distance=4):
    offending_offsets = set()
    offending_indels = set()
    for offset in mh.offsets:
        overlapping = index.query(mh.chrom, offset, distance=distance)
        if len(overlapping) > 0:
            indels = [interval.data for interval in overlapping]
            offending_indels.update(indels)
            offending_offsets.add(offset)
    if len(offending_offsets) > 0:
        max_freq = max([max(indel.alternate_frequencies) for indel in offending_indels])
        offset_data = ";".join(sorted(map(str, offending_offsets)))
        offending_indels = sorted(offending_indels, key=lambda i: (i.start, i.rsid))
        indel_data = ";".join(
            sorted(set([f"{indel.rsid}:{indel.frequency}" for indel in offending_indels]))
        )
        print(mh.name, len(mh), mh.data.Ae, max_freq, offset_data, indel_data, sep="\t")


def get_parser():
    desc = "Identify markers with allele-defining SNPs (ADSs) near common indels"
    parser = ArgumentParser(description=desc)
    parser.add_argument("markers", help="path to MicroHapDB marker definitions in CSV format")
    parser.add_argument("dbsnp", help="path to dbSNP VCF file; must be tabix-indexed (.tbi)")
    parser.add_argument(
        "--distance",
        type=int,
        default=4,
        metavar="D",
        help="flag markers with ADSs within D bp of a dbSNP indel; by default D=4",
    )
    parser.add_argument(
        "--delta",
        type=int,
        default=25,
        metavar="Δ",
        help="extend Δ bp beyond each marker when retrieving dbSNP indels; by default Δ=25",
    )
    parser.add_argument(
        "--min-freq",
        type=float,
        default=0.005,
        metavar="F",
        help="ignore indels with frequency < F; by default F=0.005 (0.5%%)",
    )
    parser.add_argument(
        "--aes",
        metavar="PATH",
        help="path to MicroHapDB Ae table in CSV format; by default, Ae values are reported as 0.0",
    )
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    markers = load_markers(args.markers, args.aes)
    main(markers, args.dbsnp, delta=args.delta, min_freq=args.min_freq, distance=args.distance)
