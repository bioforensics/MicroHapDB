# -*- coding: utf-8 -*-
#
# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse
import microhapdb
from . import lookup, marker, population, frequency
import sys
import textwrap


subparser_funcs = {
    'lookup': lookup.subparser,
    'marker': marker.subparser,
    'population': population.subparser,
    'frequency': frequency.subparser,
}

mains = {
     'lookup': lookup.main,
     'marker': marker.main,
     'population': population.main,
     'frequency': frequency.main,
}

bubbletext = r'''
≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
 __  __ _            _  _           ___  ___
|  \/  (_)__ _ _ ___| || |__ _ _ __|   \| _ )
| |\/| | / _| '_/ _ \ __ / _` | '_ \ |) | _ \
|_|  |_|_\__|_| \___/_||_\__,_| .__/___/|___/
                              |_|
≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
'''


def print_files():
    """Print the location of the installed data files on the system.

    This is helpful if the user wants to data science the TSV tables directly.
    """
    tables = (
        'marker', 'population', 'frequency', 'variantmap', 'idmap'
    )
    for table in tables:
        print(microhapdb.data_file(table + '.tsv'))


def get_parser():
    """Construct an argument parser for the command-line interface."""
    cli = argparse.ArgumentParser(
        description=bubbletext,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    cli._optionals.title = 'Configuration'
    cli.add_argument(
        '-v', '--version', action='version',
        version='MicroHapDB v{}'.format(microhapdb.__version__)
    )
    cli.add_argument(
        '-f', '--files', action='store_true', help='print data table filenames and exit'
    )

    subcommandstr = ', '.join(sorted(subparser_funcs.keys()))
    subparsers = cli.add_subparsers(dest='cmd', metavar='cmd', help=subcommandstr)
    for func in subparser_funcs.values():
        func(subparsers)

    return cli


def main(args=None):
    """MicroHapDB main method."""
    if args is None:  # pragma: no cover
        args = get_parser().parse_args()

    # If no arguments are provided, invoke --help mode
    try:
        if set([getattr(args, key) for key in vars(args)]) == set([False, None]):
            get_parser().parse_args(['-h'])
    except TypeError:
        pass

    if args.files:
        print_files()
        raise SystemExit(0)

    assert args.cmd in mains
    mainmethod = mains[args.cmd]
    mainmethod(args)
