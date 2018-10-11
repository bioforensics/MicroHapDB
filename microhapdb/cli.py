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
import pytest
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
    subparser.add_argument('query', nargs='?')

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
    subparser.add_argument('query', nargs='?')

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
    subparser.add_argument('query', nargs='?')

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
    subparser.add_argument('query', nargs='?')

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
    if args.query is None:
        result = table
    else:
        result = table.query(args.query)
    print(result.to_string())


def test_help(capsys):
    with pytest.raises(SystemExit):
        get_parser().parse_args(['-h'])
    out, err = capsys.readouterr()
    assert 'show this help message and exit' in out


def test_version(capsys):
    with pytest.raises(SystemExit):
        get_parser().parse_args(['-v'])
    out, err = capsys.readouterr()
    assert microhapdb.__version__ in out or microhapdb.__version__ in err


def test_parser():
    p = get_parser().parse_args(['locus'])
    assert p.cmd == 'locus'
    assert p.query is None

    p = get_parser().parse_args(['locus', 'Chrom == 5'])
    assert p.cmd == 'locus'
    assert p.query is 'Chrom == 5'


def test_parser_files_query(capsys):
    with pytest.raises(SystemExit) as se:
        args = get_parser().parse_args(['files', 'ID == "bogus"'])
    out, err = capsys.readouterr()
    assert 'unrecognized arguments: ID == "bogus"' in err


def test_main(capsys):
    args = get_parser().parse_args(['population'])
    main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split('\n')
    assert len(outlines) == 85, len(outlines)


def test_main_query(capsys):
    args = get_parser().parse_args(['population', 'Name.str.contains("Amer")'])
    main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split('\n')
    assert len(outlines) == 4, len(outlines)


def test_main_files(capsys):
    args = get_parser().parse_args(['files'])
    main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split('\n')
    assert len(outlines) == 4, len(outlines)
