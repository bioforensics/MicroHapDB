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
    desc = 'Retrieve marker and population records by name or identifier'
    epilog = """\
    Examples::

        microhapdb lookup rs10815466
        microhapdb lookup --format=detail mh12KK-043
    """
    epilog = dedent(epilog)
    subparser = subparsers.add_parser(
        'lookup', description=desc, epilog=epilog, formatter_class=RawDescriptionHelpFormatter,
    )
    subparser.add_argument('--format', choices=['table', 'detail'], default='table')
    subparser.add_argument('id', help='record identifier')


def main(args):
    result = microhapdb.retrieve.by_id(args.id)
    print(result.to_string(index=False))
