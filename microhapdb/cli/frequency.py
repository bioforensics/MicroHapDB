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
    subparser.add_argument('--marker', metavar='ID', help='restrict frequencies by marker')
    subparser.add_argument('--population', metavar='ID', help='restrict frequencies by population')
    subparser.add_argument('--allele', metavar='ID', help='restrict frequencies by allele')


def main(args):
    query_args = list()
    if args.marker:
        markerid = standardize_marker_ids([args.marker])
        criterion = 'Marker == "{m:s}"'.format(m=markerid.iloc[0])
        query_args.append(criterion)
    if args.population:
        popid = standardize_population_ids([args.population])
        criterion = 'Population == "{p:s}"'.format(p=popid.iloc[0])
        query_args.append(criterion)
    if args.allele:
        criterion = 'Allele == "{a:s}"'.format(a=args.allele)
        query_args.append(criterion)
    if len(query_args) > 0:
        query = ' and '.join(query_args)
        result = microhapdb.frequencies.query(query, engine='python')
    else:
        result = microhapdb.frequencies
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
