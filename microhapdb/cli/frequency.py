# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from argparse import RawDescriptionHelpFormatter
import microhapdb
from microhapdb.retrieve import standardize_marker_ids, standardize_population_ids
from numpy import float64
import pandas as pd
import sys
from textwrap import dedent


def subparser(subparsers):
    desc = microhapdb.cli.bubbletext + '\nRetrieve population allele frequencies'
    epilog = """\
    Examples::

        microhapdb frequency --marker=mh22KK-060 --population=SA000001B
        microhapdb frequency --marker=mh22KK-060 --allele=C,A
    """
    epilog = dedent(epilog)
    subparser = subparsers.add_parser(
        'frequency', description=desc, epilog=epilog, formatter_class=RawDescriptionHelpFormatter,
    )
    subparser.add_argument('--format', choices=['table', 'detail', 'mhpl8r', 'efm'], default='table')
    meg = subparser.add_mutually_exclusive_group()
    meg.add_argument('--marker', metavar='ID', nargs='+', help='restrict frequencies by marker')
    meg.add_argument('--panel', metavar='FILE', help='restrict frequencies to markers listed in FILE, one ID per line')
    subparser.add_argument('--population', metavar='ID', nargs='+', help='restrict frequencies by population')
    subparser.add_argument('--allele', metavar='ID', help='restrict frequencies by allele')


def order_haplotypes(markers, pop, frequencies):
    """Retrieve all unique haplotype designators and order by marker occurrence

    This ordering creates a pleasing arrangement of blocks of data along a single diagonal in the
    final output table, whereas other orderings result in a more chaotic organization.
    """
    haplotypes = list()
    for marker in sorted(markers):
        # print("DEBUG", frequencies.shape, len(frequencies.Marker == marker), len(frequencies.Population))
        subset = frequencies[(frequencies.Marker == marker) & (frequencies.Population == pop)]
        for haplotype in sorted(subset.Allele):
            if haplotype not in haplotypes:
                haplotypes.append(haplotype)
    print(
        f"[microhapdb] retrieved and ordered {len(haplotypes)} distinct haplotypes",
        file=sys.stderr,
    )
    return haplotypes


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


def main(args):
    query_args = list()
    result = microhapdb.frequencies
    if args.marker:
        markerids = standardize_marker_ids(args.marker)
        result = result[result.Marker.isin(markerids)]
    if args.panel:
        with open(args.panel, 'r') as fh:
            markerids = fh.read().strip().split()
            markerids = standardize_marker_ids(markerids)
        result = result[result.Marker.isin(markerids)]
    if args.population:
        popids = standardize_population_ids(args.population)
        result = result[result.Population.isin(popids)]
    if args.allele:
        result = result[result.Allele == args.allele]
    if args.format == 'table':
        print(result.to_string(index=False))
    elif args.format == 'detail':
        raise NotImplementedError('detail format not yet implemented')
    elif args.format == 'mhpl8r':
        npop = len(result.Population.unique())
        if npop > 1:
            print(f'warning: frequencies for {npop} populations recovered, expected only 1', file=sys.stderr)
        result = result[['Marker', 'Allele', 'Frequency']].rename(columns={'Allele': 'Haplotype'})
        result.to_csv(sys.stdout, sep="\t", index=False)
    elif args.format == 'efm':
        if args.population is None or len(args.population) != 1:
            raise ValueError("must specify one and only one population with --format=efm")
        result = construct_frequency_table(args.population[0], markerids)
        result.to_csv(sys.stdout)
    else:
        raise ValueError(f'unsupported view format "{args.format}"')
