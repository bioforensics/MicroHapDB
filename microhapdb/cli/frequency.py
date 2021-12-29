# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from argparse import RawDescriptionHelpFormatter
import microhapdb
from microhapdb.retrieve import standardize_marker_ids, standardize_population_ids
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
    subparser.add_argument('--format', choices=['table', 'detail', 'mhpl8r'], default='table')
    meg = subparser.add_mutually_exclusive_group()
    meg.add_argument('--marker', metavar='ID', nargs='+', help='restrict frequencies by marker')
    meg.add_argument('--panel', metavar='FILE', help='restrict frequencies to markers listed in FILE, one ID per line')
    subparser.add_argument('--population', metavar='ID', nargs='+', help='restrict frequencies by population')
    subparser.add_argument('--allele', metavar='ID', help='restrict frequencies by allele')


def main(args):
    query_args = list()
    result = microhapdb.frequencies
    if args.marker:
        markerids = standardize_marker_ids(args.marker)
        result = result[result.Marker.isin(markerids)]
    if args.panel:
        with open(args.panel, 'r') as fh:
            markerids = fh.read().strip().split()
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
    else:
        raise ValueError(f'unsupported view format "{args.format}"')
