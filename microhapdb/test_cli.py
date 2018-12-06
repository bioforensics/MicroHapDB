#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import microhapdb
from microhapdb.cli import get_parser
from microhapdb.cli import id_mode, query_mode, region_mode
import pytest


def test_main_no_args(capsys):
    with pytest.raises(SystemExit):
        args = get_parser().parse_args([])
        microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    message = 'show this help message and exit'
    assert message in out or message in err


def test_help(capsys):
    with pytest.raises(SystemExit):
        get_parser().parse_args(['-h'])
    out, err = capsys.readouterr()
    message = 'show this help message and exit'
    assert message in out or message in err


def test_version(capsys):
    with pytest.raises(SystemExit):
        get_parser().parse_args(['-v'])
    out, err = capsys.readouterr()
    assert microhapdb.__version__ in out or microhapdb.__version__ in err


def test_files(capsys):
    with pytest.raises(SystemExit):
        args = get_parser().parse_args(['--files'])
        microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    print(err)
    outlines = out.strip().split('\n')
    assert len(outlines) == 6


def test_parser():
    cli = get_parser().parse_args(['--table', 'locus'])
    assert cli.table == 'locus'
    assert cli.query is None
    assert cli.region is None

    cli = get_parser().parse_args(['--table', 'locus', '--region', 'chr5'])
    assert cli.table == 'locus'
    assert cli.query is None
    assert cli.region == 'chr5'

    cli = get_parser().parse_args(
        ['--table', 'population', '--query', 'Name.str.contains("Amer")']
    )
    assert cli.table == 'population'
    assert cli.query == 'Name.str.contains("Amer")'
    assert cli.region is None


def test_main(capsys):
    args = get_parser().parse_args(['--table', 'population'])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split('\n')
    assert len(outlines) == 96 + 1  # 96 populations + 1 header line


def test_main_region_mode(capsys):
    arglist = ['--region', 'chr15:52192400-52192500', '--table', 'variant']
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split('\n')
    print(out)
    assert len(outlines) == 5 + 1  # 5 variants + 1 header line


def test_main_id_mode(capsys):
    arglist = ['--id', 'rs10815466']
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    line = '34996  MHDBV000034997    GRCh38  chr9    680713     A,G  dbSNP151'
    assert line in out


def test_main_warnings(capsys):
    arglist = [
        '--query', 'Name.str.contains("Afr")', '--table', 'population',
        '--id', 'rs10815466',
    ]
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    assert 'ignoring "id" parameter in "query" mode' in err

    arglist = [
        '--query', 'Name.str.contains(",")', '--table', 'population',
        '--region', 'chr5:10000000-20000000',
    ]
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    assert 'ignoring "region" parameter in "query" mode' in err

    arglist = ['--id', 'rs10815466', '--table', 'population']
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    assert 'ignoring "table" parameter in "id" mode' in err
