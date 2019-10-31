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
from microhapdb.util import data_file
import pytest
from tempfile import NamedTemporaryFile


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


def test_main_pop_noargs(capsys):
    args = get_parser().parse_args(['population'])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split('\n')
    assert len(outlines) == 1 + 96 + 3 + 1 + 1 + 1  # 1 header line + 96 ALFRED + 3 LOVD + 1 Linköping + 1 NRIPS + 1 Chen2019


def test_main_pop_detail(capsys):
    args = get_parser().parse_args(['population', '--format=detail', 'Koreans'])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    assert '876 total allele frequencies available' in out


def test_main_pop_query(capsys):
    args = get_parser().parse_args(['population', '--query', 'Source == "10.1016/j.fsigen.2018.05.008"'])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    testout = '''
               ID    Name                        Source
 MHDBP-3dab7bdd14  Africa  10.1016/j.fsigen.2018.05.008
 MHDBP-936bc36f79    Asia  10.1016/j.fsigen.2018.05.008
 MHDBP-383d86606a      NL  10.1016/j.fsigen.2018.05.008
'''
    print(out)
    assert testout.strip() == out.strip()


def test_main_marker_noargs(capsys):
    args = get_parser().parse_args(['marker'])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split('\n')
    assert len(outlines) == 198 + 15 + 40 + 26 + 1 + 11


def test_main_marker_detail(capsys):
    args = get_parser().parse_args(['marker', '--format=detail', 'mh01CP-008'])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    assert '>mh01CP-008\nGACATCACGCCACTGCT\n' in out


def test_main_marker_query(capsys):
    args = get_parser().parse_args(['marker', '--query', 'Chrom == "chr19"'])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    testout = '''
       Name          PermID Reference  Chrom                                       Offsets   AvgAe                          Source
 mh19KK-056  MHDBM-d6ff8635    GRCh38  chr19                               4852124,4852324  2.4143                          ALFRED
 mh19CP-007  MHDBM-49dbcc57    GRCh38  chr19                    14310739,14310772,14310780  3.0813                          ALFRED
  mh19NH-23  MHDBM-dd72537b    GRCh38  chr19                    22052723,22052774,22052817     NaN  10.1016/j.legalmed.2015.06.003
 mh19KK-299  MHDBM-8cbeb11c    GRCh38  chr19  22546697,22546748,22546779,22546810,22546850  3.8989                          ALFRED
  mh19AT-47  MHDBM-8f439540    GRCh38  chr19                    22546697,22546748,22546779  1.4537                   ISFG2019:P597
 mh19KK-301  MHDBM-2069446a    GRCh38  chr19           50938487,50938502,50938526,50938550  1.8143                          ALFRED
 mh19KK-057  MHDBM-eb558c37    GRCh38  chr19                    51654948,51655025,51655062  2.1923                          ALFRED
'''
    assert testout.strip() == out.strip()


def test_main_marker_fasta(capsys):
    args = get_parser().parse_args(['marker', '--format=fasta', 'mh13CP-010', 'mh08PK-46625'])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    testout = '''
>mh08PK-46625 PermID=MHDBM-840756f3 GRCh38:chr8:1194352-1194371 variants=115,119,127,134
TGCTGGCAAGTTGGAAACACAGGCTCTGCGATTTTGAGAGTGAACCTGCAAGAGACAAGCAGACGTTGGCAGTGCCGCGT
CCCGGCTGGTGGAGGGAGCCCGGATGCCTGGCAGACAGTCAGTGGTCGGTTGGCGGCCGGCCCACATAAGGGCACCATGC
TCACCGTGTCTAGGCAGAGCTGGAGGCTCCTCCTGCCCAGGGCGGCCTCCAGGTGGGGAGGACGGCAGAGCTTCCCTCAG
TCCCACTTTC
>mh13CP-010 PermID=MHDBM-13233c9a GRCh38:chr13:29218044-29218076 variants=109,120,141
AATAAGACCTGGTCTCCACAAAGAAATTTTAAAAATTAGCTGGGCTTGGTGATGCATGCCTGTAGTCCCAGCTACTGAGG
CTGAGGCAGGAGTATTCCTTGAGTCCAGGAGGTCATGGCTGCAGTGAGTTATGATTGTGCCGTCATACTCCAGCCTGAAC
AAAAGAGTGAGACCTTGTCCCTCCCCGCCAAAACCAAACCAAAACAAAACAAAACAAAAAAAAAACACCTAAAAACCCCA
GTGTTTACAGT
'''
    print(out)
    assert testout.strip() == out.strip()


def test_main_marker_region_mode(capsys):
    arglist = ['marker', '--region', 'chr15']
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split('\n')
    print(out)
    assert len(outlines) == 13 + 1  # 13 markers + 1 header line


def test_main_marker_region_mode_failure(capsys):
    arglist = ['marker', '--region', 'chr15:']
    args = get_parser().parse_args(arglist)
    with pytest.raises(ValueError, match=r'cannot parse region "chr15:"'):
        microhapdb.cli.main(args)


def test_main_marker_panel(capsys):
    with NamedTemporaryFile() as panelfile:
        with open(panelfile.name, 'w') as fh:
            for marker in ['mh05KK-058', 'mh06KK-101', 'mh20KK-035']:
                print(marker, file=fh)
        arglist = ['marker', '--panel', panelfile.name]
        args = get_parser().parse_args(arglist)
        microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    testout = '''
       Name          PermID Reference  Chrom                     Offsets   AvgAe  Source
 mh06KK-101  MHDBM-8a2c760e    GRCh38   chr6         170280714,170280900  1.7068  ALFRED
 mh05KK-058  MHDBM-d6c594d2    GRCh38  chr15  28120284,28120471,28120586  2.1016  ALFRED
 mh20KK-035  MHDBM-92f3685a    GRCh38  chr20             2088698,2088728  2.0937  ALFRED
'''
    print(terminal.out)
    assert testout.strip() == terminal.out.strip()


def test_main_marker_panel_query_conflict(capsys):
    with NamedTemporaryFile() as panelfile:
        with open(panelfile.name, 'w') as fh:
            for marker in ['mh05KK-058', 'mh06KK-101', 'mh20KK-035']:
                print(marker, file=fh)
        arglist = ['marker', '--panel', panelfile.name, '--query', 'PermID == "MHDBM-8a2c760e"']
        args = get_parser().parse_args(arglist)
        microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    assert 'ignoring user-supplied marker IDs in --query mode' in terminal.err


def test_main_marker_panel_region_conflict(capsys):
    with NamedTemporaryFile() as panelfile:
        with open(panelfile.name, 'w') as fh:
            for marker in ['mh05KK-058', 'mh06KK-101', 'mh20KK-035']:
                print(marker, file=fh)
        arglist = ['marker', '--panel', panelfile.name, '--region=chr12']
        args = get_parser().parse_args(arglist)
        microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    assert 'ignoring user-supplied marker IDs in --region mode' in terminal.err


@pytest.mark.parametrize('pop,marker,allele,numrows', [
    ('--population=Swedish', None, None, 138),
    ('--population=SA000009J', '--marker=mh13KK-218', None, 15),
    (None, '--marker=mh13KK-218', '--allele=C,T,C,T', 97),
    (None, '--marker=mh14PK-72639', None, 46),
    (None, None, None, 82807)
])
def test_main_frequency_by_pop(pop, marker, allele, numrows, capsys):
    testargs = (pop, marker, allele)
    arglist = ['frequency'] + [arg for arg in testargs if arg is not None]
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    outlines = terminal.out.strip().split('\n')
    assert len(outlines) == numrows


@pytest.mark.parametrize('panel', [
    'alpha',
    'beta',
])
def test_main_panel(panel, capsys):
    arglist = ['marker', '--panel', panel, '--format=fasta']
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    testout = data_file('tests/panel-' + panel + '.fasta')
    print(terminal.out)
    with open(testout, 'r') as fh:
        assert fh.read().strip() == terminal.out.strip()


def test_lookup(capsys):
    arglist = ['lookup', 'rs10815466']
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    assert 'mh09KK-033  MHDBM-8458b727    GRCh38  chr9         680713,680762,680790  2.8101         ALFRED' in out
    assert ' mh09AT-15  MHDBM-b46abf2e    GRCh38  chr9  680713,680762,680767,680790  2.8991  ISFG2019:P597' in out
