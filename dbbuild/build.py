#!/usr/bin/env python
#
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
from lib import SourceIndex
from pathlib import Path
import sys


def main(source_path, dbsnp_path, chain_path, exclusions=["Auton2015"], check_only=False):
    validate_paths(dbsnp_path, chain_path)
    if check_only:
        return
    index = SourceIndex(source_path, dbsnp_path, chain_path, exclude=exclusions)
    index.interval_check()
    index.update_marker_names()
    index.markers.to_csv("marker.csv", index=False)
    index.indels.to_csv("indels.csv", index=False)
    index.frequencies.to_csv(
        "frequency.csv.gz", index=False, float_format="%.5f", compression="gzip"
    )
    index.populations.to_csv("population.csv", index=False)
    index.merges.to_csv("merged.csv", index=False)
    print(index)


def validate_paths(dbsnp_path, chain_path):
    paths = list()
    for version in (37, 38):
        for extension in ("vcf.gz", "vcf.gz.tbi", "rsidx"):
            path = Path(dbsnp_path) / f"dbSNP_GRCh{version}.{extension}"
            paths.append(path)
    paths.append(Path(chain_path) / "hg19ToHg38.over.chain.gz")
    paths.append(Path(chain_path) / "hg38ToHg19.over.chain.gz")
    files_present = [p.is_file() for p in paths]
    print("-" * 60, "[Auxiliary data file check]\n", "Present  Path", sep="\n", file=sys.stderr)
    for path, present in zip(paths, files_present):
        print(f"{present!s:8} {path}", file=sys.stderr)
    print("-" * 60, file=sys.stderr)
    if sum(files_present) < len(paths):
        missing = [str(path) for path, present in zip(paths, files_present) if present is False]
        raise FileNotFoundError(",".join(missing))


def get_parser():
    parser = ArgumentParser(description="MicroHapDB database build procedure")
    parser.add_argument("dbsnp_path")
    parser.add_argument("chain_path")
    parser.add_argument(
        "--sources",
        default="sources",
        metavar="PATH",
        help="path to directory containing primary data sources; by default PATH=sources/",
    )
    parser.add_argument(
        "--exclude",
        nargs="+",
        metavar="SRC",
        default=["Auton2015"],
        help="one or more sources to exclude from the database build; by default SRC=auton2015",
    )
    parser.add_argument(
        "--check", action="store_true", help="perform auxiliary data file check only and exit"
    )
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    main(
        args.sources,
        args.dbsnp_path,
        args.chain_path,
        exclusions=args.exclude,
        check_only=args.check,
    )
