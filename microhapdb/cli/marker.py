# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from argparse import RawDescriptionHelpFormatter
import microhapdb
from microhapdb.marker import print_detail, print_table
from textwrap import dedent


def subparser(subparsers):
    desc = microhapdb.cli.bubbletext + '\nRetrieve marker records by identifier or query'
    epilog = """\
    Examples::

        microhapdb marker mh01NK-001
        microhapdb marker --format=detail --minlen=225 MHDBM-dc55cd9e
        microhapdb marker --region=chr18:1-25000000
        microhapdb marker --query='Source == "ALFRED"'
        microhapdb marker --query='Name.str.contains("PK")'
    """
    epilog = dedent(epilog)
    subparser = subparsers.add_parser(
        'marker', description=desc, epilog=epilog, formatter_class=RawDescriptionHelpFormatter,
    )
    subparser.add_argument('--format', choices=['table', 'detail'], default='table')
    subparser.add_argument(
        '--delta', metavar='D', type=int, default=25, help='extend D nucleotides beyond the '
        'marker extent when computing amplicon boundaries (detail format only); by default D=25'
    )
    subparser.add_argument(
        '--min-length', metavar='L', type=int, default=250, help='minimum amplicon length (detail '
        'format only); by default L=250'
    )
    subparser.add_argument(
        '-r', '--region', metavar='RGN', help='restrict results to the '
        'specified genomic region; format chrX:YYYY-ZZZZZ'
    )
    subparser.add_argument(
        '--query', metavar='STRING', help='Retrieve records using a Pandas-style query'
    )
    subparser.add_argument('id', nargs='*', help='marker identifier')


def main(args):
    if args.query:
        result = microhapdb.markers.query(args.query)
    elif args.region:
        result = microhapdb.retrieve.by_region(args.region)
    elif len(args.id) > 0:
        idents = microhapdb.retrieve.standardize_marker_ids(args.id)
        result = microhapdb.markers[microhapdb.markers.Name.isin(idents)]
    else:
        result = microhapdb.markers
    view = print_table if args.format == 'table' else print_detail
    view(result, delta=args.delta, minlen=args.min_length)
