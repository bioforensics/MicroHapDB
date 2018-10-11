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
import textwrap


tables = {
    'allelefreq': 'allelefreqs',
    'locus': 'loci',
    'population': 'populations',
    'variant': 'variants',
    'files': None,
}


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
    subcommandstr = '", "'.join(sorted(list(tables.keys())))
    parser = argparse.ArgumentParser(
        description=bubbletext,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser._optionals.title = 'Global arguments'
    parser.add_argument(
        '-v', '--version', action='version',
        version='MicroHapDB v{}'.format(microhapdb.__version__)
    )
    subparsers = parser.add_subparsers(dest='cmd', metavar='cmd',
                                       help='"' + subcommandstr + '"')

    desc = 'List or retrieve allele frequencies for 148 loci across 84 populations.'
    desc = textwrap.dedent(desc)
    epilog = """\
    Examples::

        microhapdb allelefreq
        microhapdb allelefreq 'Locus == "SI664601X"'
        microhapdb allelefreq 'Locus == "SI664601X" & Allele == "C,T,G"'
        microhapdb allelefreq 'Locus == "SI664601X" & Allele == "C,T,G" & Population == "SA000007H"'
    """
    epilog = textwrap.dedent(epilog)
    subparser = subparsers.add_parser(
        'allelefreq', description=desc, epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparser.add_argument('query', nargs='*')

    desc = 'List or retrieve data on microhaplotype loci'
    epilog = """\
    Examples::

        microhapdb locus
        microhapdb locus 'ID == "SI664572E"'
        microhapdb locus 'Name == "mh05KK-122"'
        microhapdb locus 'Variants.str.contains("rs10408594")'
        microhapdb locus 'Chrom == 2 & Start > 200000000 & End < 300000000'
    """
    epilog = textwrap.dedent(epilog)
    subparser = subparsers.add_parser(
        'locus', description=desc, epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparser.add_argument('query', nargs='*')

    desc = 'List or retrieve data on the 84 populations included in MicroHapDB.'
    epilog = """\
    Examples::

        microhapdb population
        microhapdb population 'ID == "SA000936S"'
        microhapdb population 'Name.str.contains("Afr")'
    """
    epilog = textwrap.dedent(epilog)
    subparser = subparsers.add_parser(
        'population', description=desc, epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparser.add_argument('query', nargs='*')

    desc = 'List or retrieve data on microhaplotype variants.'
    epilog = """\
    Examples::

        microhapdb variant
        microhapdb variant 'ID == "rs10815466"'
        microhapdb variant 'Chrom == 7'
    """
    epilog = textwrap.dedent(epilog)
    subparser = subparsers.add_parser(
        'variant', description=desc, epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparser.add_argument('query', nargs='*')

    desc = 'List MicroHapDB data files'
    subparser = subparsers.add_parser('files', description=desc)

    return parser


def main(args=None):
    if args is None:  # pragma: no cover
        args = get_parser().parse_args()

    if args.cmd is None:  # pragma: no cover
        get_parser().parse_args(['-h'])

    assert args.cmd in tables
    if args.cmd == 'files':
        for datatype in ('allele', 'locus', 'population', 'variant'):
            print(microhapdb.data_file(datatype + '.tsv'))
        return
    table = getattr(microhapdb, tables[args.cmd])
    if len(args.query) == 0:
        result = table
    elif len(args.query) == 1:
        result = table.query(args.query[0])
    else:
        raise ValueError('multiple queries unsupported')
    print(result.to_string())
