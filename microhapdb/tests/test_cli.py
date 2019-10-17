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
import pytest


def test_main_no_args(capsys):
    with pytest.raises(SystemExit):
        args = get_parser().parse_args([])
        microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    message = 'show this help message and exit'
    assert message in terminal.out or message in terminal.err


def test_help(capsys):
    with pytest.raises(SystemExit):
        get_parser().parse_args(['-h'])
    terminal = capsys.readouterr()
    message = 'show this help message and exit'
    assert message in terminal.out or message in terminal.err


def test_version(capsys):
    with pytest.raises(SystemExit):
        get_parser().parse_args(['-v'])
    terminal = capsys.readouterr()
    assert microhapdb.__version__ in terminal.out or microhapdb.__version__ in terminal.err


def test_files(capsys):
    with pytest.raises(SystemExit):
        args = get_parser().parse_args(['--files'])
        microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    print(err)
    outlines = out.strip().split('\n')
    assert len(outlines) == 7


def test_parser():
    cli = get_parser().parse_args(['marker'])
    assert cli.cmd == 'marker'
    assert cli.query is None
    assert cli.region is None
    assert cli.id == []

    cli = get_parser().parse_args(['marker', '--region', 'chr5'])
    assert cli.cmd == 'marker'
    assert cli.query is None
    assert cli.region == 'chr5'
    assert cli.id == []

    cli = get_parser().parse_args(['marker', 'mh04CP-003', 'mh17KK-014'])
    assert cli.cmd == 'marker'
    assert cli.query is None
    assert cli.region is None
    assert cli.id == ['mh04CP-003', 'mh17KK-014']

    cli = get_parser().parse_args(['population', '--query', 'Name.str.contains("Amer")'])
    assert cli.cmd == 'population'
    assert cli.query == 'Name.str.contains("Amer")'
    assert cli.id == []


def test_main(capsys):
    args = get_parser().parse_args(['population'])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split('\n')
    assert len(outlines) == 1 + 96 + 3 + 1  # 1 header line + 96 ALFRED + 3 LOVD + 1 Linköping


def test_main_region_mode(capsys):
    arglist = ['marker', '--region', 'chr15']
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split('\n')
    print(out)
    assert len(outlines) == 13 + 1  # 13 markers + 1 header line


def test_main_region_mode_failure(capsys):
    arglist = ['marker', '--region', 'chr15:']
    args = get_parser().parse_args(arglist)
    with pytest.raises(ValueError, match=r'cannot parse region "chr15:"'):
        microhapdb.cli.main(args)


def test_lookup(capsys):
    arglist = ['lookup', 'rs10815466']
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    assert 'mh09KK-033  MHDBM-8458b727    GRCh38  chr9         680713,680762,680790  2.8101         ALFRED' in out
    assert ' mh09AT-15  MHDBM-b46abf2e    GRCh38  chr9  680713,680762,680767,680790  2.8991  ISFG2019:P597' in out
