# -------------------------------------------------------------------------------------------------
# Copyright (c) 2019, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from argparse import RawDescriptionHelpFormatter
import microhapdb
from microhapdb import Marker, Population
from numpy import float64
import pandas as pd
import sys
from textwrap import dedent
from warnings import warn


def main(args):
    markers = args.marker
    if args.panel:
        with open(args.panel, "r") as fh:
            markers = fh.read().strip().split()
    result = apply_filters(markers, args.population, args.allele)
    display(result, args.format, args.population)


def apply_filters(markers=None, populations=None, allele=None):
    result = microhapdb.frequencies
    if markers:
        markers = Marker.standardize_ids(markers)
        result = result[result.Marker.isin(markers)]
    if populations:
        populations = Population.standardize_ids(populations)
        result = result[result.Population.isin(populations)]
    if allele:
        result = result[result.Allele == allele]
    return result


def display(result, view_format, population):
    if view_format == "table":
        result.to_csv(sys.stdout, sep="\t", index=False)
    elif view_format == "detail":
        raise NotImplementedError("detail format not yet implemented")
    elif view_format == "mhpl8r":
        npop = len(result.Population.unique())
        if npop > 1:
            warn(f"frequencies for {npop} populations recovered, expected only 1", UserWarning)
        result = result[["Marker", "Allele", "Frequency"]].rename(columns={"Allele": "Haplotype"})
        result.to_csv(sys.stdout, sep="\t", index=False)
    elif view_format == "efm":
        if population is None or len(population) != 1:
            raise ValueError("must specify one and only one population with --format=efm")
        result = construct_frequency_table(args.population[0], markerids)
        result.to_csv(sys.stdout)
    else:
        raise ValueError(f'unsupported view format "{view_format}"')


def subparser(subparsers):
    desc = microhapdb.cli.bubbletext + "\nRetrieve population allele frequencies"
    epilog = """\
    Examples::

        microhapdb frequency --marker=mh22KK-060 --population=SA000001B
        microhapdb frequency --marker=mh22KK-060 --allele=C,A
    """
    epilog = dedent(epilog)
    subparser = subparsers.add_parser(
        "frequency",
        description=desc,
        epilog=epilog,
        formatter_class=RawDescriptionHelpFormatter,
    )
    subparser.add_argument(
        "--format", choices=["table", "detail", "mhpl8r", "efm"], default="table"
    )
    meg = subparser.add_mutually_exclusive_group()
    meg.add_argument("--marker", metavar="ID", nargs="+", help="restrict frequencies by marker")
    meg.add_argument(
        "--panel",
        metavar="FILE",
        help="restrict frequencies to markers listed in FILE, one ID per line",
    )
    subparser.add_argument(
        "--population", metavar="ID", nargs="+", help="restrict frequencies by population"
    )
    subparser.add_argument("--allele", metavar="ID", help="restrict frequencies by allele")


def construct_frequency_table(pop, panel):
    frequencies = microhapdb.frequencies
    frequencies = frequencies[frequencies.Marker.isin(panel)]
    marker_indices = {marker: i for i, marker in enumerate(panel)}
    haplotypes = order_haplotypes(panel, pop, frequencies)
    table = pd.DataFrame(index=haplotypes, columns=sorted(panel), dtype=float64)
    table.index.name = "Allele"
    for haplotype in haplotypes:
        hapfreqs = [None] * len(panel)
        subset = frequencies[(frequencies.Allele == haplotype) & (frequencies.Population == pop)]
        for j, row in subset.iterrows():
            index = marker_indices[row.Marker]
            hapfreqs[index] = row.Frequency
        table.loc[haplotype] = hapfreqs
    nrows, ncols = table.shape
    print(
        f"[microhapdb] constructed frequency table for {nrows} haplotypes and {ncols} markers",
        file=sys.stderr,
    )
    return table.round(3)


def order_haplotypes(markers, pop, frequencies):
    """Retrieve all unique haplotype designators and order by marker occurrence

    This ordering creates a pleasing arrangement of blocks of data along a single diagonal in the
    final output table, whereas other orderings result in a more chaotic organization.
    """
    haplotypes = list()
    for marker in sorted(markers):
        subset = frequencies[(frequencies.Marker == marker) & (frequencies.Population == pop)]
        for haplotype in sorted(subset.Allele):
            if haplotype not in haplotypes:
                haplotypes.append(haplotype)
    print(
        f"[microhapdb] retrieved and ordered {len(haplotypes)} distinct haplotypes",
        file=sys.stderr,
    )
    return haplotypes
