#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse
import microhapdb
from microhapdb.retrieve import query_mode, id_mode, region_mode
import sys
import textwrap


def print_files():
    tables = (
        'allele', 'locus', 'population', 'variant', 'variantmap', 'idmap'
    )
    for table in tables:
        print(microhapdb.data_file(table + '.tsv'))


def get_parser():
    bubbletext = r'''
≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
 __  __ _            _  _           ___  ___
|  \/  (_)__ _ _ ___| || |__ _ _ __|   \| _ )
| |\/| | / _| '_/ _ \ __ / _` | '_ \ |) | _ \
|_|  |_|_\__|_| \___/_||_\__,_| .__/___/|___/
                              |_|
≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
'''
    epilog = """\
    Examples::

        microhapdb --id mh12KK-043
        microhapdb --id rs10815466
        microhapdb --id SA000936S
        microhapdb --table population
        microhapdb --table locus --region chr7
        microhapdb --table variant --region chr19:4000000-5000000
        microhapdb --table allele --query 'Locus == "MHDBL000047" and Allele == "C,T,G" and Population == "MHDBP000006"'
    """
    epilog = textwrap.dedent(epilog)
    cli = argparse.ArgumentParser(
        description=bubbletext, epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    cli._optionals.title = 'Configuration'
    cli.add_argument(
        '-v', '--version', action='version',
        version='MicroHapDB v{}'.format(microhapdb.__version__)
    )
    cli.add_argument(
        '-f', '--files', action='store_true', help='print data table '
        'filenames and exit'
    )
    cli.add_argument(
        '-t', '--table', choices=['variant', 'locus', 'population', 'allele'],
        metavar='TBL', help='restrict results to the specified data table; '
        'must be one of "variant", "locus", "population", or "allele"'
    )
    cli.add_argument(
        '-r', '--region', metavar='RGN', help='restrict results to the '
        'specified genomic region; format chrX:YYYY-ZZZZZ'
    )
    cli.add_argument(
        '--id', metavar='ID', help='query data tables using a dbSNP ID, an '
        'ALFRED ID/name, or an internal MicroHapDB ID'
    )
    cli.add_argument(
        '-q', '--query', metavar='QRY', help='Invoke a Pandas-style query; '
        'must specify table to query with the `-t|--table` flag'
    )

    return cli


def main(args=None):
    if args is None:  # pragma: no cover
        args = get_parser().parse_args()

    if set([getattr(args, key) for key in vars(args)]) == set([False, None]):
        get_parser().parse_args(['-h'])

    elif args.files:
        print_files()

    elif args.query:
        for attr in ('region', 'id'):
            if getattr(args, attr):
                msg = 'ignoring "{}" parameter in "query" mode'.format(attr)
                print('[MicroHapDB] WARNING:', msg, file=sys.stderr)
        query_mode(args.table, args.query)

    elif args.id:
        if args.table:
            message = 'ignoring "table" parameter in "id" mode'
            print('[MicroHapDB] WARNING:', message, file=sys.stderr)
        id_mode(args.id)

    elif args.region:
        region_mode(args.region, args.table)

    else:
        print(microhapdb.tables[args.table].to_string())
