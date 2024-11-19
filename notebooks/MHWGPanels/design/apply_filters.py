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
from collections import Counter
import matplotlib
from matplotlib import pyplot as plt
import pandas as pd
import sys
from upsetplot import UpSet


def main(markers, masks, pass_stream=sys.stdout, fail_stream=None):
    filter_markers(masks, fail_stream=fail_stream)
    filter_accounting(markers, masks)
    passed_markers = markers[~markers.Name.isin(masks.filtered_markers)]
    passed_markers.to_csv(pass_stream, index=False)


def filter_markers(masks, fail_stream=None):
    fail_modes = Counter()
    if fail_stream:
        print(FilteredMarker.header(), file=fail_stream)
    for filtered_marker in masks:
        if fail_stream:
            print(filtered_marker, file=fail_stream)
        fail_modes[filtered_marker.fail_mode] += 1
    print(*fail_modes.most_common(), sep="\n", file=sys.stderr)


def filter_accounting(markers, masks):
    fail_modes = Counter()
    filtered_markers = set()
    for filtered_marker in masks:
        filtered_markers.add(filtered_marker)
        fail_modes[filtered_marker.filtering] += 1
    prepassed_markers = markers[~markers.Name.isin(masks.filtered_markers)]
    prepassed_length_fail = prepassed_markers[prepassed_markers.Extent > 260]
    fail_modes[(False, True, False, False, False)] = len(prepassed_length_fail)
    fail_tuples = sorted(fail_modes.items())
    keys = (k for k, v in fail_tuples)
    counts = (v for k, v in fail_tuples)
    index = pd.MultiIndex.from_tuples(
        keys,
        names=[
            "Indel Interference",
            "Length",
            "Low Complexity",
            "Repetitive Content",
            "STR Proximity",
        ],
    )
    s = pd.Series(counts, index=index)
    backend = matplotlib.get_backend()
    plt.switch_backend("Agg")
    fig = plt.figure(figsize=(12, 6), dpi=300)
    UpSet(s, show_counts=True).plot(fig=fig)
    plt.savefig("upset.png", bbox_inches="tight")
    plt.switch_backend(backend)


class MaskSet(dict):
    def __init__(self, indel_path, lowcomplex_path, repeat_path, str_path, maxlen=260):
        self.masks = {
            "indel": pd.read_csv(indel_path, sep="\t").set_index("Marker"),
            "lowcomplex": pd.read_csv(lowcomplex_path, sep="\t").set_index("Marker"),
            "repeat": pd.read_csv(repeat_path, sep="\t").set_index("Marker"),
            "str": pd.read_csv(str_path, sep="\t").set_index("Marker"),
        }
        self.maxlen = maxlen

    def __iter__(self):
        for marker in sorted(self.filtered_markers):
            yield self.query(marker)

    @property
    def filtered_markers(self):
        filtered = set()
        for mask in self.masks.values():
            filtered.update(mask.index)
        return filtered

    def query(self, marker_name):
        results = {}
        for key, mask in self.masks.items():
            results[key] = mask.loc[marker_name] if marker_name in mask.index else None
        return FilteredMarker(results, maxlen=self.maxlen)


class FilteredMarker:
    def __init__(self, mask_data, maxlen=260):
        self.maxlen = maxlen
        self.mask_data = mask_data
        for key, entry in mask_data.items():
            if entry is None:
                continue
            self.marker = entry.name
            self.extent = entry.Extent
            self.ae = entry.Ae
            break
        else:
            raise ValueError("empty mask data?")

    def __len__(self):
        return self.extent

    @property
    def fail_mode(self):
        failures = set()
        for key, mask in self.mask_data.items():
            if mask is not None:
                failures.add(key)
        if len(self) > self.maxlen:
            failures.add("length")
        return "".join(f"[{key}]" for key in sorted(failures))

    @staticmethod
    def header():
        return "\t".join(
            (
                "Marker",
                "Extent",
                "Ae",
                "FailMode",
                "MaxAltFreq",
                "LowComplexRepeats",
                "Repeats",
                "STRs",
            )
        )

    def __str__(self):
        return f"{self.marker}\t{self.extent}\t{self.ae}\t{self.fail_mode}\t{self.max_alt_freq}\t{self.low_complex_repeats}\t{self.repeats}\t{self.strs}"

    @property
    def max_alt_freq(self):
        if self.mask_data["indel"] is None:
            return ""
        return self.mask_data["indel"].MaxAltFreq

    @property
    def low_complex_repeats(self):
        if self.mask_data["lowcomplex"] is None:
            return ""
        return self.mask_data["lowcomplex"].Repeats

    @property
    def repeats(self):
        if self.mask_data["repeat"] is None:
            return ""
        return self.mask_data["repeat"].Repeats

    @property
    def strs(self):
        if self.mask_data["str"] is None:
            return ""
        return self.mask_data["str"].Locus

    @property
    def filtering(self):
        return (
            self.mask_data["indel"] is not None,
            len(self) > 260,
            self.mask_data["lowcomplex"] is not None,
            self.mask_data["repeat"] is not None,
            self.mask_data["str"] is not None,
        )


def get_parser():
    parser = ArgumentParser(description="Discard markers based on the provided mask(s)")
    parser.add_argument("markers", help="path to MicroHapDB marker definitions in CSV format")
    parser.add_argument("indel_mask", help="path to CSV file with indel mask")
    parser.add_argument("lowcomplex_mask", help="path to CSV file with low-complexity mask")
    parser.add_argument("repeat_mask", help="path to CSV file with repeat mask")
    parser.add_argument("str_mask", help="path to CSV file with forensic STR mask")
    parser.add_argument(
        "--maxlen",
        type=int,
        metavar="L",
        default=260,
        help="filter out markers longer than L bp; by default L=260",
    )
    parser.add_argument(
        "--passed",
        metavar="FILE",
        help="output file for markers passing filters; default is terminal (stdout)",
    )
    parser.add_argument("--failed", metavar="FILE", help="output file for markers failing filters")
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    markers = pd.read_csv(args.markers)
    masks = MaskSet(args.indel_mask, args.lowcomplex_mask, args.repeat_mask, args.str_mask)
    pass_stream = sys.stdout if args.passed is None else open(args.passed, "w")
    fail_stream = None if args.failed is None else open(args.failed, "w")
    main(markers, masks, pass_stream=pass_stream, fail_stream=fail_stream)
    if args.passed:
        pass_stream.close()
    if args.failed:
        fail_stream.close()
