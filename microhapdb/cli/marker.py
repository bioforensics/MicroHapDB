# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from argparse import RawDescriptionHelpFormatter
import microhapdb
from microhapdb.marker import print_detail, print_table, print_fasta
from textwrap import dedent
import sys


def subparser(subparsers):
    desc = microhapdb.cli.bubbletext + '\nRetrieve marker records by identifier or query'
    epilog = """\
    Examples::

        microhapdb marker mh01NK-001
        microhapdb marker --format=fasta mh13KK-218 mh04CP-002 mh02AT-05
        microhapdb marker --format=fasta --panel mypanel.txt
        microhapdb marker --format=detail --min-length=125 MHDBM-dc55cd9e
        microhapdb marker --region=chr18:1-25000000
        microhapdb marker --query='Source == "ALFRED"'
        microhapdb marker --query='Name.str.contains("PK")'
    """
    epilog = dedent(epilog)
    subparser = subparsers.add_parser(
        'marker', description=desc, epilog=epilog, formatter_class=RawDescriptionHelpFormatter,
    )
    subparser.add_argument('--format', choices=['table', 'detail', 'fasta'], default='table')
    subparser.add_argument(
        '--notrunc', dest='trunc', action='store_false', default=True,
        help='disable truncation of tabular results'
    )
    subparser.add_argument(
        '--delta', metavar='D', type=int, default=25, help='extend D nucleotides beyond the '
        'marker extent when computing amplicon boundaries (detail and fasta format only); by '
        'default D=25'
    )
    subparser.add_argument(
        '--min-length', metavar='L', type=int, default=250, help='minimum amplicon length (detail '
        'and fasta format only); by default L=250'
    )
    subparser.add_argument(
        '-p', '--panel', metavar='FILE', help='file containing a list of marker names/identifiers,'
        ' one per line'
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
        if len(args.id) > 0 or args.panel is not None:
            warning = 'WARNING: ignoring user-supplied marker IDs in --query mode'
            print('[MicroHapDB::marker]', warning, file=sys.stderr)
        result = microhapdb.markers.query(args.query)
    elif args.region:
        if len(args.id) > 0 or args.panel is not None:
            warning = 'WARNING: ignoring user-supplied marker IDs in --region mode'
            print('[MicroHapDB::marker]', warning, file=sys.stderr)
        result = microhapdb.retrieve.by_region(args.region)
    elif len(args.id) > 0 or args.panel is not None:
        idents = args.id
        if args.panel:
            if hasattr(microhapdb.panel, args.panel):
                func = getattr(microhapdb.panel, args.panel)
                idents.extend(func())
            else:
                with open(args.panel, 'r') as fh:
                    for line in fh:
                        idents.append(line.strip())
        idents = microhapdb.retrieve.standardize_marker_ids(idents)
        result = microhapdb.markers[microhapdb.markers.Name.isin(idents)]
    else:
        result = microhapdb.markers
    viewfuncs = {
        'table': print_table,
        'detail': print_detail,
        'fasta': print_fasta
    }
    view = viewfuncs[args.format]
    view(result, delta=args.delta, minlen=args.min_length, trunc=args.trunc)
