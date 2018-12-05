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
from sys import stderr
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

    if args.files:
        print_files()
        return

    if set([getattr(args, key) for key in vars(args)]) == set([False, None]):
        get_parser().parse_args(['-h'])

    if args.query:
        for attr in ('region', 'id'):
            if getattr(args, attr):
                msg = 'ignoring "{}" parameter in "query" mode'.format(attr)
                print('[MicroHapDB] WARNING:', message, file=stderr)
        query_mode(args.table, args.query)

    elif args.id:
        if args.table:
            message = 'ignoring "table" parameter in "id" mode'
            print('[MicroHapDB] WARNING:', message, file=stderr)
        id_mode(args.id)

    elif args.region:
        region_mode(args.region, args.table)

    else:
        print(microhapdb.tables[args.table].to_string())



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
    assert len(outlines) == 84, len(outlines)


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
