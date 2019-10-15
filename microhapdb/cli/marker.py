#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from argparse import RawDescriptionHelpFormatter
import microhapdb
from textwrap import dedent


def subparser(subparsers):
    desc = 'Retrieve marker records by identifier or query'
    epilog = """\
    Examples::

        microhapdb marker mh01NK-001
        microhapdb marker --format=detail MHDBM-dc55cd9e
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
    print(result.to_string(index=False))
