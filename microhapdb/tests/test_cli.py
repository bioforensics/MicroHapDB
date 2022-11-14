# -------------------------------------------------------------------------------------------------
# Copyright (c) 2018, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from io import StringIO
import microhapdb
from microhapdb.cli import get_parser
import pandas
from pkg_resources import resource_filename
import pytest
from tempfile import NamedTemporaryFile



def data_file(path):
    return resource_filename("microhapdb", f"tests/data/{path}")


def test_main_no_args(capsys):
    with pytest.raises(SystemExit):
        args = get_parser().parse_args([])
        microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    message = "show this help message and exit"
    assert message in terminal.out or message in terminal.err


def test_help(capsys):
    with pytest.raises(SystemExit):
        get_parser().parse_args(["-h"])
    terminal = capsys.readouterr()
    message = "show this help message and exit"
    assert message in terminal.out or message in terminal.err


def test_version(capsys):
    with pytest.raises(SystemExit):
        get_parser().parse_args(["-v"])
    terminal = capsys.readouterr()
    assert microhapdb.__version__ in terminal.out or microhapdb.__version__ in terminal.err


def test_files(capsys):
    with pytest.raises(SystemExit):
        args = get_parser().parse_args(["--files"])
        microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split("\n")
    assert len(outlines) == 7


def test_parser():
    cli = get_parser().parse_args(["marker"])
    assert cli.cmd == "marker"
    assert cli.query is None
    assert cli.region is None
    assert cli.id == []

    cli = get_parser().parse_args(["marker", "--region", "chr5"])
    assert cli.cmd == "marker"
    assert cli.query is None
    assert cli.region == "chr5"
    assert cli.id == []

    cli = get_parser().parse_args(["marker", "mh04CP-003", "mh17KK-014"])
    assert cli.cmd == "marker"
    assert cli.query is None
    assert cli.region is None
    assert cli.id == ["mh04CP-003", "mh17KK-014"]

    cli = get_parser().parse_args(["population", "--query", 'Name.str.contains("Amer")'])
    assert cli.cmd == "population"
    assert cli.query == 'Name.str.contains("Amer")'
    assert cli.id == []


def test_main_pop_noargs(capsys):
    args = get_parser().parse_args(["population"])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split("\n")
    assert len(outlines) == 1 + 109


def test_main_pop_detail(capsys):
    args = get_parser().parse_args(["population", "--format=detail", "Koreans"])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    assert "876 total allele frequencies available" in out


def test_main_pop_query(capsys):
    args = get_parser().parse_args(
        ["population", "--query", 'Source == "10.1016/j.fsigen.2018.05.008"']
    )
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    testout = """
              ID   Name                       Source
MHDBP-3dab7bdd14 Africa 10.1016/j.fsigen.2018.05.008
MHDBP-936bc36f79   Asia 10.1016/j.fsigen.2018.05.008
MHDBP-383d86606a     NL 10.1016/j.fsigen.2018.05.008
"""
    print(out)
    assert testout.strip() == out.strip()


def test_main_marker_noargs(capsys):
    args = get_parser().parse_args(["marker"])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split("\n")
    assert len(outlines) == 1 + 634


def test_main_marker_notrunc(capsys):
    args = get_parser().parse_args(["marker", "--notrunc"])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split("\n")
    assert len(outlines) == 1 + 634


def test_main_marker_detail(capsys):
    args = get_parser().parse_args(["marker", "--format=detail", "mh01CP-008"])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    assert ">mh01CP-008\nGACATCACGCCACTGCT\n" in out


def test_main_marker_query(capsys):
    args = get_parser().parse_args(["marker", "--query", 'Chrom == "chr19"'])
    microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    exp_out = """
        Name         PermID Reference Chrom                                                                          Offsets     Ae     In     Fst                         Source
mh19USC-19pA MHDBM-2d713eab    GRCh38 chr19                                                             561778,561795,561815 2.7453 0.0870  0.0733   10.1016/j.fsigen.2019.102213
  mh19KK-056 MHDBM-d6ff8635    GRCh38 chr19                                                                  4852124,4852324 2.4391 0.0773  0.1760                         ALFRED
 mh19SHY-001 MHDBM-bc5d0364    GRCh38 chr19                  7698912,7698928,7698991,7699013,7699030,7699043,7699090,7699096 4.8141 0.2843  0.0985     10.1007/s00414-020-02483-x
  mh19CP-007 MHDBM-49dbcc57    GRCh38 chr19                                                       14310739,14310772,14310780 3.0813 0.0466  0.0776                         ALFRED
mh19USC-19pB MHDBM-76427274    GRCh38 chr19                                     16040864,16040894,16040899,16040921,16040929 3.5107 0.1647  0.0731   10.1016/j.fsigen.2019.102213
 mh19ZHA-006 MHDBM-5f597339    GRCh38 chr19                            20579862,20579863,20579879,20579892,20579916,20579924 3.8739 0.1681 -0.0032   10.1016/j.fsigen.2020.102255
   mh19NH-23 MHDBM-dd72537b    GRCh38 chr19                                                       22052723,22052774,22052817 1.9380 0.0414 -0.0124 10.1016/j.legalmed.2015.06.003
mh19KKCS-299 MHDBM-a70896aa    GRCh38 chr19                   22546697,22546702,22546748,22546779,22546810,22546829,22546850    NaN    NaN     NaN   10.1016/j.fsigen.2020.102275
  mh19KK-299 MHDBM-8cbeb11c    GRCh38 chr19                                     22546697,22546748,22546779,22546810,22546850 4.1592 0.2335  0.1102                         ALFRED
   mh19AT-47 MHDBM-8f439540    GRCh38 chr19                                                       22546697,22546748,22546779 2.4025 0.1298  0.1170                  ISFG2019:P597
 mh19ZHA-007 MHDBM-afb5b4d1    GRCh38 chr19                                              28397315,28397408,28397447,28397455 3.7262 0.1551 -0.0007            10.1098/rsos.191937
mh19USC-19qA MHDBM-f757e745    GRCh38 chr19                                              33273771,33273785,33273811,33273816 3.3219 0.1127  0.0880   10.1016/j.fsigen.2019.102213
  mh19KK-301 MHDBM-2069446a    GRCh38 chr19                                              50938487,50938502,50938526,50938550 1.9707 0.2673  0.1698                         ALFRED
mh19KKCS-300 MHDBM-bc8b7213    GRCh38 chr19                   50947786,50947789,50947790,50947830,50947876,50947877,50947967    NaN    NaN     NaN   10.1016/j.fsigen.2020.102275
  mh19KK-057 MHDBM-eb558c37    GRCh38 chr19                                                       51654948,51655025,51655062 2.3266 0.0667  0.0428                         ALFRED
 mh19ZHA-009 MHDBM-346aa500    GRCh38 chr19                                     53129072,53129075,53129133,53129203,53129249 3.8739 0.1405 -0.0241            10.1098/rsos.191937
mh19USC-19qB MHDBM-7b40359b    GRCh38 chr19                                                       53714387,53714389,53714413 4.0756 0.1368  0.0163   10.1016/j.fsigen.2019.102213
 mh19SHY-002 MHDBM-2b448141    GRCh38 chr19 55588420,55588451,55588465,55588474,55588492,55588517,55588545,55588551,55588584 3.2215 0.2184  0.0049     10.1007/s00414-020-02483-x
"""
    obs_out = terminal.out
    print(obs_out)
    assert exp_out.strip() == obs_out.strip()


def test_main_marker_fasta_default_delta(capsys):
    args = get_parser().parse_args(["marker", "--format=fasta", "mh01CP-016", "mh06PK-24844"])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    testout = """
>mh01CP-016 PermID=MHDBM-021e569a GRCh38:chr1:55559012-55559057 variants=18,56,62 Xref=SI664876L
TGGCACACAACAAGTGCTTATAATGAAAGCATTAGTGAGTAAAAGAGTGATCCCTGGCTTTGAACTCCCTCTAAGTGTAC
C
>mh06PK-24844 PermID=MHDBM-aa39cbba GRCh38:chr6:13861392-13861447 variants=13,20,35,42,51,55,59,60,61,67
AGGAAGAAAGTGATTACATCCAAACGTGAGCAGGAGGAAACTCGGAACATACTGTTTTTAAGAACTAGTATCACTAGAGT
T
"""
    print(out)
    assert testout.strip() == out.strip()


def test_main_marker_fasta_long_delta(capsys):
    args = get_parser().parse_args(
        [
            "marker",
            "--format=fasta",
            "--delta=25",
            "--min-length=250",
            "mh13CP-010",
            "mh08PK-46625",
        ]
    )
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    testout = """
>mh08PK-46625 PermID=MHDBM-840756f3 GRCh38:chr8:1194352-1194372 variants=115,119,127,134
TGCTGGCAAGTTGGAAACACAGGCTCTGCGATTTTGAGAGTGAACCTGCAAGAGACAAGCAGACGTTGGCAGTGCCGCGT
CCCGGCTGGTGGAGGGAGCCCGGATGCCTGGCAGACAGTCAGTGGTCGGTTGGCGGCCGGCCCACATAAGGGCACCATGC
TCACCGTGTCTAGGCAGAGCTGGAGGCTCCTCCTGCCCAGGGCGGCCTCCAGGTGGGGAGGACGGCAGAGCTTCCCTCAG
TCCCACTTTC
>mh13CP-010 PermID=MHDBM-13233c9a GRCh38:chr13:29218044-29218077 variants=109,120,141
AATAAGACCTGGTCTCCACAAAGAAATTTTAAAAATTAGCTGGGCTTGGTGATGCATGCCTGTAGTCCCAGCTACTGAGG
CTGAGGCAGGAGTATTCCTTGAGTCCAGGAGGTCATGGCTGCAGTGAGTTATGATTGTGCCGTCATACTCCAGCCTGAAC
AAAAGAGTGAGACCTTGTCCCTCCCCGCCAAAACCAAACCAAAACAAAACAAAACAAAAAAAAAACACCTAAAAACCCCA
GTGTTTACAGT
"""
    print(out)
    assert testout.strip() == out.strip()


def test_main_marker_region_mode(capsys):
    arglist = ["marker", "--region", "chr15"]
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split("\n")
    print(out)
    assert len(outlines) == 24 + 1  # markers + 1 header line


def test_main_marker_region_mode_failure(capsys):
    arglist = ["marker", "--region", "chr15:"]
    args = get_parser().parse_args(arglist)
    with pytest.raises(ValueError, match=r'cannot parse region "chr15:"'):
        microhapdb.cli.main(args)


def test_main_marker_panel(capsys):
    with NamedTemporaryFile() as panelfile:
        with open(panelfile.name, "w") as fh:
            for marker in ["mh15KK-058", "mh06KK-101", "mh20KK-035"]:
                print(marker, file=fh)
        arglist = ["marker", "--panel", panelfile.name]
        args = get_parser().parse_args(arglist)
        microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    testout = """
      Name         PermID Reference Chrom                    Offsets     Ae     In    Fst Source
mh06KK-101 MHDBM-8a2c760e    GRCh38  chr6        170280714,170280900 1.6705 0.2296 0.2590 ALFRED
mh15KK-058 MHDBM-d6c594d2    GRCh38 chr15 28120284,28120471,28120586 2.2110 0.3799 0.0862 ALFRED
mh20KK-035 MHDBM-92f3685a    GRCh38 chr20            2088698,2088728 2.1328 0.2104 0.2180 ALFRED
"""
    print(terminal.out)
    assert testout.strip() == terminal.out.strip()


def test_main_marker_panel_plus_query(capsys):
    with NamedTemporaryFile() as panelfile:
        with open(panelfile.name, "w") as fh:
            for marker in ["mh05KK-058", "mh06KK-101", "mh20KK-035"]:
                print(marker, file=fh)
        arglist = ["marker", "--panel", panelfile.name, "--query", 'PermID == "MHDBM-8a2c760e"']
        args = get_parser().parse_args(arglist)
        microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    outlines = terminal.out.strip().split("\n")
    assert len(outlines) == 2


@pytest.mark.parametrize(
    "mode,offsets1,offsets2",
    [
        ("symmetric", "variants=70,92,104", "variants=47,82,128"),
        ("5", "variants=120,142,154", "variants=73,108,154"),
        ("3", "variants=20,42,54", "variants=20,55,101"),
    ],
)
def test_main_marker_extendmode(mode, offsets1, offsets2, capsys):
    arglist = [
        "marker",
        "--format",
        "fasta",
        "--extend-mode",
        mode,
        "--delta",
        "20",
        "--min-length",
        "175",
        "mh01USC-1pD",
        "mh22NH-27",
    ]
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    assert offsets1 in terminal.out
    assert offsets2 in terminal.out


@pytest.mark.parametrize("mode", ["4", "6", "NotARealMode"])
def test_main_marker_extendmode_bad(mode, capsys):
    arglist = [
        "marker",
        "--format",
        "fasta",
        "--extend-mode",
        mode,
        "--delta",
        "20",
        "--min-length",
        "175",
        "mh01USC-1pD",
        "mh22NH-27",
    ]
    with pytest.raises(SystemExit):
        args = get_parser().parse_args(arglist)
    terminal = capsys.readouterr()
    assert "invalid str_to_extend_mode value" in terminal.err


def test_main_marker_view_bad():
    arglist = [
        "marker",
        "--format",
        "fasta",
        "--extend-mode",
        "5",
        "--delta",
        "20",
        "--min-length",
        "175",
        "mh01USC-1pD",
        "mh22NH-27",
    ]
    args = get_parser().parse_args(arglist)
    args.format = "html"
    with pytest.raises(ValueError, match=r'unsupported view format "html"'):
        microhapdb.cli.main(args)


@pytest.mark.parametrize(
    "pop,marker,allele,numrows",
    [
        ("--population=Swedish", None, None, 138),
        ("--population=SA000009J", "--marker=mh13KK-218", None, 15),
        (None, "--marker=mh13KK-218", "--allele=C,T,C,T", 97),
        (None, "--marker=mh14PK-72639", None, 217),
        (None, None, None, 149248),
    ],
)
def test_main_frequency_by_pop(pop, marker, allele, numrows, capsys):
    testargs = (pop, marker, allele)
    arglist = ["frequency"] + [arg for arg in testargs if arg is not None]
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    outlines = terminal.out.strip().split("\n")
    assert len(outlines) == numrows


@pytest.mark.parametrize(
    "panel",
    [
        "alpha",
        "beta",
    ],
)
def test_main_panel(panel, capsys):
    arglist = [
        "marker",
        "--panel",
        panel + "_legacy",
        "--format=fasta",
        "--delta=25",
        "--min-length=250",
    ]
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    observed = terminal.out.strip()
    testout = data_file(f"panel-{panel}.fasta")
    with open(testout, "r") as fh:
        expected = fh.read().strip()
    assert observed == expected


def test_lookup(capsys):
    arglist = ["lookup", "rs10815466"]
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    assert (
        "mh09KK-033 MHDBM-8458b727    GRCh38  chr9        680713,680762,680790 2.9343 0.1008 0.0504        ALFRED"
        in out
    )
    assert (
        " mh09AT-15 MHDBM-b46abf2e    GRCh38  chr9 680713,680762,680767,680790 2.9471 0.1160 0.0602 ISFG2019:P597"
        in out
    )


def test_ae_pop(capsys):
    arglist = ["marker", "--region=chr18:1-25000000", "--ae-pop=CDX"]
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    exp_out = """
        Name         PermID Reference Chrom                                                                                 Offsets     Ae     In    Fst                       Source
 mh18SHY-001 MHDBM-b8a01c31    GRCh38 chr18 1952476,1952501,1952508,1952522,1952549,1952565,1952604,1952642,1952646,1952655,1952656 3.8490 0.2426 0.0230   10.1007/s00414-020-02483-x
mh18PK-87558 MHDBM-1e5374f1    GRCh38 chr18                                         1960542,1960557,1960561,1960566,1960582,1960588 1.9767 0.1325 0.0494 10.1016/j.fsigen.2018.05.008
mh18USC-18pA MHDBM-56dfa93b    GRCh38 chr18                                                 5280017,5280020,5280070,5280071,5280093 1.8402 0.2130 0.1818 10.1016/j.fsigen.2019.102213
  mh18CP-005 MHDBM-a85754d3    GRCh38 chr18                                                         8892864,8892893,8892896,8892907 3.3873 0.0904 0.0059                       ALFRED
 mh18ZBF-002 MHDBM-3702ab96    GRCh38 chr18                                                     11900702,11900723,11900734,11900778 3.2413 0.2019 0.0864       10.1002/elps.201900451
 mh18ZHA-004 MHDBM-db649dba    GRCh38 chr18                                                     14315931,14315952,14315958,14316046 3.1131 0.1949 0.0118 10.1016/j.fsigen.2020.102255
  mh18KK-285 MHDBM-ea520d26    GRCh38 chr18                                                     24557354,24557431,24557447,24557489 2.6188 0.1721 0.0836                       ALFRED
   mh18AT-38 MHDBM-db09ec41    GRCh38 chr18                                                              24557431,24557447,24557489 2.6144 0.1419 0.0837                ISFG2019:P597
"""
    obs_out = terminal.out
    print(obs_out)
    assert exp_out.strip() == obs_out.strip()


def test_ae_pop_bad_pop():
    arglist = ["marker", "--ae-pop=ABC", "mh18USC-18pA"]
    args = get_parser().parse_args(arglist)
    with pytest.raises(ValueError, match=r'no Ae data for population "ABC"'):
        microhapdb.cli.main(args)


def test_hg37(capsys):
    arglist = ["marker", "--region=chr18:1-25000000", "--GRCh37"]
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    exp_out = """
        Name         PermID Reference Chrom                                                                                 Offsets     Ae     In    Fst                       Source
 mh18SHY-001 MHDBM-b8a01c31    GRCh37 chr18 1952477,1952502,1952509,1952523,1952550,1952566,1952605,1952643,1952647,1952656,1952657 4.4183 0.2426 0.0230   10.1007/s00414-020-02483-x
mh18PK-87558 MHDBM-1e5374f1    GRCh37 chr18                                         1960543,1960558,1960562,1960567,1960583,1960589 2.3659 0.1325 0.0494 10.1016/j.fsigen.2018.05.008
mh18USC-18pA MHDBM-56dfa93b    GRCh37 chr18                                                 5280016,5280019,5280069,5280070,5280092 3.4330 0.2130 0.1818 10.1016/j.fsigen.2019.102213
  mh18CP-005 MHDBM-a85754d3    GRCh37 chr18                                                         8892862,8892891,8892894,8892905 3.6722 0.0904 0.0059                       ALFRED
 mh18ZBF-002 MHDBM-3702ab96    GRCh37 chr18                                                     11900701,11900722,11900733,11900777 3.2962 0.2019 0.0864       10.1002/elps.201900451
 mh18ZHA-004 MHDBM-db649dba    GRCh37 chr18                                                     14315930,14315951,14315957,14316045 2.9579 0.1949 0.0118 10.1016/j.fsigen.2020.102255
  mh18KK-285 MHDBM-ea520d26    GRCh37 chr18                                                     22137318,22137395,22137411,22137453 2.7524 0.1721 0.0836                       ALFRED
   mh18AT-38 MHDBM-db09ec41    GRCh37 chr18                                                              22137395,22137411,22137453 2.7093 0.1419 0.0837                ISFG2019:P597
  mh18CP-003 MHDBM-6fdf83f9    GRCh37 chr18                                                              23068395,23068425,23068433 3.1124 0.1061 0.0183                       ALFRED
 mh18ZBF-001 MHDBM-742f43e2    GRCh37 chr18                                                                       24223721,24223752 2.9546 0.1796 0.0047       10.1002/elps.201900451
"""
    obs_out = terminal.out
    print(obs_out)
    assert exp_out.strip() == obs_out.strip()


def test_hg37_detail(capsys):
    arglist = ["marker", "--format=detail", "--GRCh37", "mh17USC-17pA"]
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    observed = terminal.out
    expected = """
--------------------------------------------------------------[ MicroHapDB ]----
mh17USC-17pA    a.k.a MHDBM-d5646523

Marker Definition (GRCh37)
    Marker extent
        - chr17:3821918-3821988 (70 bp)
    Target locus
        - chr17:3918614-3918704 (90 bp)
    Constituent variants
        - chromosome offsets: 3821918,3821952,3821987
        - marker offsets: 0,34,69
        - target offsets: 10,44,79
        - cross-references: rs4995288, rs4995289, rs9904113
    Observed haplotypes
        - C,C,A
        - C,C,C
        - C,T,A
        - C,T,C
        - T,C,A
        - T,C,C
        - T,T,A
        - T,T,C


--[ Core Marker Sequence ]--
>mh17USC-17pA
CGTCTCATTTGGGGATCTTATCATATCCACAGTGTACCCCAGGGACCTACATTTATTTTCCTGCTGCTGA


--[ Marker Target Sequence with MH alleles (haplotypes) ]--
          *                                 *                                  *
CCTGGAGCACCGTCTCATTTGGGGATCTTATCATATCCACAGTGTACCCCAGGGACCTACATTTATTTTCCTGCTGCTGATAACTCCACA
..........C.................................C..................................A..........
..........C.................................C..................................C..........
..........C.................................T..................................A..........
..........C.................................T..................................C..........
..........T.................................C..................................A..........
..........T.................................C..................................C..........
..........T.................................T..................................A..........
..........T.................................T..................................C..........
--------------------------------------------------------------------------------
"""
    assert observed.strip() == expected.strip()


def test_hg37_ae_pop(capsys):
    arglist = ["marker", "--region=chr18:50000000-80000000", "--GRCh37", "--ae-pop=PJL"]
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    exp_out = """
        Name         PermID Reference Chrom                                                                                            Offsets     Ae     In     Fst                       Source
mh18USC-18qB MHDBM-14fcada5    GRCh37 chr18                                                                         50547498,50547528,50547540 3.8565 0.1449  0.0500 10.1016/j.fsigen.2019.102213
mh18USC-18qC MHDBM-6bf74efc    GRCh37 chr18                                                                63842523,63842541,63842557,63842562 3.5570 0.1436 -0.0049 10.1016/j.fsigen.2019.102213
 mh18SHY-002 MHDBM-721f5fbc    GRCh37 chr18 73681402,73681414,73681416,73681430,73681439,73681497,73681525,73681529,73681563,73681574,73681588 4.8436 0.2109  0.0714   10.1007/s00414-020-02483-x
mh18KKCS-293 MHDBM-350bd971    GRCh37 chr18                                     76089731,76089843,76089884,76089885,76089906,76089944,76089967    NaN    NaN     NaN 10.1016/j.fsigen.2020.102275
  mh18KK-293 MHDBM-13ed6da8    GRCh37 chr18                                                                76089885,76089906,76089944,76089967 2.6445 0.2495  0.0837                       ALFRED
   mh18AT-39 MHDBM-13ed6da8    GRCh37 chr18                                                                76089885,76089906,76089944,76089967 2.6445 0.2495  0.0837                ISFG2019:P597
"""
    obs_out = terminal.out
    print(obs_out)
    assert exp_out.strip() == obs_out.strip()


def test_marker_offsets_cli(capsys):
    arglist = [
        "marker",
        "--format=offsets",
        "--delta=25",
        "--min-length=200",
        "mh03AT-09",
        "mh11KK-180",
        "mh13KK-217",
        "mh07USC-7qC",
    ]
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    result = pandas.read_csv(StringIO(terminal.out), sep="\t")
    assert result.shape == (15, 4)
    observed = list(result.Offset)
    expected = [85, 114, 66, 95, 122, 123, 134, 25, 145, 203, 218, 25, 65, 179, 217]
    assert observed == expected
    observed = list(result.OffsetHg38)
    expected = [
        131927127,
        131927156,
        151821663,
        151821692,
        151821719,
        151821720,
        151821731,
        1669560,
        1669680,
        1669738,
        1669753,
        46291794,
        46291834,
        46291948,
        46291986,
    ]
    assert observed == expected


def test_marker_offsets_37_cli(capsys):
    arglist = [
        "marker",
        "--GRCh37",
        "--format=offsets",
        "mh01USC-1qA",
        "mh02USC-2qA",
        "mh03USC-3qA",
    ]
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    result = pandas.read_csv(StringIO(terminal.out), sep="\t")
    assert result.shape == (10, 4)
    assert list(result.columns) == ["Marker", "Offset", "Chrom", "OffsetHg37"]
    assert list(result.OffsetHg37)[:5] == [167126967, 167126986, 103092502, 103092512, 103092574]


def test_mhpl8r(capsys):
    arglist = ["frequency", "--marker", "mh02USC-2pA", "--population", "JPT", "--format", "mhpl8r"]
    args = microhapdb.cli.get_parser().parse_args(arglist)
    microhapdb.cli.frequency.main(args)
    terminal = capsys.readouterr()
    result = pandas.read_csv(StringIO(terminal.out), sep="\t")
    assert result.shape == (4, 3)
    assert result.Haplotype.iloc[0] == "A,A,G,A"
    assert result.Frequency.iloc[0] == pytest.approx(0.005)


def test_mhpl8r_panel(capsys):
    with NamedTemporaryFile(mode="wt") as tempfile:
        print("mh02USC-2pA", "mh08USC-8qA", "mh17USC-17qA", sep="\n", file=tempfile, flush=True)
        arglist = [
            "frequency",
            "--panel",
            tempfile.name,
            "--population",
            "GBR",
            "--format",
            "mhpl8r",
        ]
        args = microhapdb.cli.get_parser().parse_args(arglist)
        microhapdb.cli.frequency.main(args)
    terminal = capsys.readouterr()
    result = pandas.read_csv(StringIO(terminal.out), sep="\t")
    print(result)
    assert result.shape == (13, 3)
    assert result.Haplotype.iloc[7] == "A,G,T"
    assert result.Frequency.iloc[7] == pytest.approx(0.429)


def test_mhpl8r_multi_pop():
    arglist = [
        "frequency",
        "--marker",
        "mh02USC-2pA",
        "mh08USC-8qA",
        "mh17USC-17qA",
        "--population",
        "CLM",
        "GIH",
        "ASW",
        "CEU",
        "--format",
        "mhpl8r",
    ]
    args = microhapdb.cli.get_parser().parse_args(arglist)
    message = r"frequencies for 4 populations recovered, expected only 1"
    with pytest.warns(UserWarning, match=message):
        microhapdb.cli.frequency.main(args)


def test_efm(capsys):
    arglist = [
        "frequency",
        "--population=CEU",
        "--format=efm",
        "--marker",
        "mh01USC-1pD",
        "mh17USC-17pA",
        "mh15USC-15qA",
    ]
    args = microhapdb.cli.get_parser().parse_args(arglist)
    microhapdb.cli.frequency.main(args)
    terminal = capsys.readouterr()
    result = pandas.read_csv(StringIO(terminal.out))
    print(result)
    assert result.shape == (13, 4)
    assert result["Allele"].iloc[3] == "C,T,C"
    assert result["mh01USC-1pD"].iloc[3] == pytest.approx(0.101)
    assert pandas.isna(result["mh15USC-15qA"].iloc[3])
    assert result["mh17USC-17pA"].iloc[3] == pytest.approx(0.02)


def test_efm_multi_pop():
    arglist = [
        "frequency",
        "--population",
        "CEU",
        "IBS",
        "--format=efm",
        "--marker",
        "mh01USC-1pD",
        "mh17USC-17pA",
        "mh15USC-15qA",
    ]
    args = microhapdb.cli.get_parser().parse_args(arglist)
    with pytest.raises(
        ValueError, match="must specify one and only one population with --format=efm"
    ):
        microhapdb.cli.frequency.main(args)


def test_bad_format():
    arglist = ["frequency", "--marker", "mh02USC-2pA", "--population", "JPT", "--format", "detail"]
    args = microhapdb.cli.get_parser().parse_args(arglist)
    with pytest.raises(NotImplementedError):
        microhapdb.cli.frequency.main(args)
    args.format = "BoGuS"
    with pytest.raises(ValueError, match=r'unsupported view format "BoGuS"'):
        microhapdb.cli.frequency.main(args)
