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
from repeats import main as flag_repeats
import sys


def main(
    source_path, dbsnp_path, chain_path, rmsk_path, exclusions=["Auton2015"], check_only=False
):
    validate_paths(dbsnp_path, chain_path, rmsk_path)
    if check_only:
        return
    index = SourceIndex(source_path, dbsnp_path, chain_path, exclude=exclusions)
    index.interval_check()
    index.update_marker_names()
    index.markers.to_csv("marker.csv", index=False)
    index.indels.to_csv("indels.csv", index=False)
    frequencies = index.frequencies
    frequencies = cleanup_frequencies(frequencies)
    frequencies.to_csv("frequency.csv.gz", index=False, float_format="%.5f", compression="gzip")
    index.populations.to_csv("population.csv", index=False)
    index.merges.to_csv("merged.csv", index=False)
    repeats = flag_repeats(Path(rmsk_path) / "rmsk.txt.gz", "marker.csv", delta=25)
    repeats.to_csv("repeats.csv", index=False)
    print(index)


def validate_paths(dbsnp_path, rmsk_path, chain_path):
    paths = list()
    for version in (37, 38):
        for extension in ("vcf.gz", "vcf.gz.tbi", "rsidx"):
            path = Path(dbsnp_path) / f"dbSNP_GRCh{version}.{extension}"
            paths.append(path)
    paths.append(Path(dbsnp_path) / "refsnp-merged.csv.gz")
    paths.append(Path(chain_path) / "hg19ToHg38.over.chain.gz")
    paths.append(Path(chain_path) / "hg38ToHg19.over.chain.gz")
    paths.append(Path(rmsk_path) / "rmsk.txt.gz")
    files_present = [p.is_file() for p in paths]
    print("-" * 60, "[Auxiliary data file check]\n", "Present  Path", sep="\n", file=sys.stderr)
    for path, present in zip(paths, files_present):
        print(f"{present!s:8} {path}", file=sys.stderr)
    print("-" * 60, file=sys.stderr)
    if sum(files_present) < len(paths):
        missing = [str(path) for path, present in zip(paths, files_present) if present is False]
        raise FileNotFoundError(",".join(missing))


def cleanup_frequencies(freq):
    freq["NumVars"] = freq.Allele.apply(lambda x: x.count("|") + 1)
    freq.loc[(freq.Marker == "mh01NK-001") & (freq.Source == "Kidd2018"), "Marker"] = "mh01NH-01.v2"
    freq.loc[(freq.Marker == "mh01NK-001") & (freq.Source == "Turchi2019"), "Marker"] = "mh01NH-01.v2"
    freq.loc[(freq.Marker == "mh01NK-001") & (freq.Source == "Gandotra2020"), "Marker"] = "mh01NH-01.v2"
    freq.loc[(freq.Marker == "mh01NK-001") & (freq.Source == "Staadig2021"), "Marker"] = "mh01NH-01.v1"
    freq.loc[(freq.Marker.str.startswith(("mh05KK-023", "mh05KK-020"))) & (freq.Source == "Kidd2018") & (freq.NumVars == 3), "Marker"] = "mh05KK-023.v1"
    freq.loc[(freq.Marker.str.startswith(("mh05KK-023", "mh05KK-020"))) & (freq.Source == "Kidd2018") & (freq.NumVars == 4), "Marker"] = "mh05KK-023.v2"
    freq.loc[(freq.Marker.str.startswith(("mh05KK-023", "mh05KK-020"))) & (freq.Source == "Turchi2019"), "Marker"] = "mh05KK-023.v1"
    freq.loc[(freq.Marker.str.startswith(("mh05KK-023", "mh05KK-020"))) & (freq.Source == "Gandotra2020"), "Marker"] = "mh05KK-023.v3"
    freq.loc[(freq.Marker.str.startswith("mh05KK-120")) & (freq.Source == "Kidd2018") & (freq.NumVars == 3), "Marker"] = "mh05KK-120.v1"
    freq.loc[(freq.Marker.str.startswith("mh05KK-120")) & (freq.Source == "Kidd2018") & (freq.NumVars == 4), "Marker"] = "mh05KK-120.v2"
    freq = freq.drop(columns="NumVars")
    freq["Count"] = freq["Count"].astype("Int16")
    return freq


def get_parser():
    parser = ArgumentParser(description="MicroHapDB database build procedure")
    parser.add_argument("dbsnp_path")
    parser.add_argument("chain_path")
    parser.add_argument("rmsk_path")
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
        args.rmsk_path,
        exclusions=args.exclude,
        check_only=args.check,
    )
