# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from argparse import RawDescriptionHelpFormatter
import microhapdb
from microhapdb.retrieve import standardize_marker_ids, standardize_population_ids
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
    subparser.add_argument('--format', choices=['table', 'detail'], default='table')
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
    print(result.to_string(index=False))
