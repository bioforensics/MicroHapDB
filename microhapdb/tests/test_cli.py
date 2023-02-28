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
    assert len(outlines) == 4


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
    assert len(outlines) == 1 + 99


def test_main_pop_detail(capsys):
    args = get_parser().parse_args(["population", "--format=detail", "Koreans"])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    assert "876 total allele frequencies available" in out


def test_main_pop_query(capsys):
    args = get_parser().parse_args(["population", "--query", 'Source == "vanderGaag2018"'])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    testout = """
              ID   Name         Source
MHDBP-3dab7bdd14 Africa vanderGaag2018
MHDBP-936bc36f79   Asia vanderGaag2018
MHDBP-383d86606a     NL vanderGaag2018
"""
    print(out)
    assert testout.strip() == out.strip()


def test_main_marker_noargs(capsys):
    args = get_parser().parse_args(["marker"])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split("\n")
    assert len(outlines) == 1 + 753


def test_main_marker_notrunc(capsys):
    args = get_parser().parse_args(["marker", "--notrunc"])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split("\n")
    assert len(outlines) == 1 + 753


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
         Name  NumVars  Extent Chrom    Start      End    Ae                   Source
 mh19USC-19pA        3      38 chr19   561779   561816 2.913           delaPuente2020
   mh19KK-056        2     201 chr19  4852125  4852325 2.594                 Kidd2018
  mh19SHY-001        8     185 chr19  7698913  7699097 6.347                   Wu2021
   mh19CP-007        3      42 chr19 14310740 14310781 3.254                 Kidd2018
 mh19USC-19pB        5      66 chr19 16040865 16040930 3.799           delaPuente2020
  mh19ZHA-006        6      63 chr19 20579863 20579925 3.167                  Sun2020
    mh19NH-23        3      95 chr19 22052724 22052818 2.024              Hiroaki2015
mh19KK-299.v1        5     154 chr19 22546698 22546851 4.060      Kidd2018;Turchi2019
mh19KK-299.v2        7     154 chr19 22546698 22546851 4.073             Gandotra2020
mh19KK-299.v4       10     182 chr19 22546698 22546879 4.073              Pakstis2021
mh19KK-299.v3        3      63 chr19 22546749 22546811 3.603              Staadig2021
  mh19ZHA-007        4     141 chr19 28397316 28397456 4.428      Kureshi2020;Sun2020
 mh19USC-19qA        4      46 chr19 33273772 33273817 3.523           delaPuente2020
   mh19KK-301        4      64 chr19 50938488 50938551 2.624      Kidd2018;Turchi2019
   mh19KK-300        7     182 chr19 50947787 50947968 5.821 Gandotra2020;Pakstis2021
   mh19KK-057        3     115 chr19 51654949 51655063 2.539      Kidd2018;Turchi2019
  mh19ZHA-009        5     178 chr19 53129073 53129250 4.347      Kureshi2020;Sun2020
 mh19USC-19qB        3      27 chr19 53714388 53714414 4.933           delaPuente2020
  mh19SHY-002        9     165 chr19 55588421 55588585 3.613                   Wu2021
"""
    obs_out = terminal.out
    print(obs_out)
    assert exp_out.strip() == obs_out.strip()


def test_main_marker_fasta_default_delta(capsys):
    args = get_parser().parse_args(["marker", "--format=fasta", "mh01CP-016", "mh06PK-24844"])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    expected = """
>mh01CP-016 GRCh38:chr1:55558994-55559075 variants=18,56,62
TGGCACACAACAAGTGCTTATAATGAAAGCATTAGTGAGTAAAAGAGTGATCCCTGGCTTTGAACTCCCTCTAAGTGTAC
C
>mh06PK-24844 GRCh38:chr6:13861379-13861460 variants=13,20,35,42,51,55,59,60,61,67
AGGAAGAAAGTGATTACATCCAAACGTGAGCAGGAGGAAACTCGGAACATACTGTTTTTAAGAACTAGTATCACTAGAGT
T
"""
    observed = out
    assert expected.strip() == observed.strip()


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
    observed = """
>mh08PK-46625 GRCh38:chr8:1194237-1194487 variants=115,119,127,134
TGCTGGCAAGTTGGAAACACAGGCTCTGCGATTTTGAGAGTGAACCTGCAAGAGACAAGCAGACGTTGGCAGTGCCGCGT
CCCGGCTGGTGGAGGGAGCCCGGATGCCTGGCAGACAGTCAGTGGTCGGTTGGCGGCCGGCCCACATAAGGGCACCATGC
TCACCGTGTCTAGGCAGAGCTGGAGGCTCCTCCTGCCCAGGGCGGCCTCCAGGTGGGGAGGACGGCAGAGCTTCCCTCAG
TCCCACTTTC
>mh13CP-010 GRCh38:chr13:29217935-29218186 variants=109,120,141
AATAAGACCTGGTCTCCACAAAGAAATTTTAAAAATTAGCTGGGCTTGGTGATGCATGCCTGTAGTCCCAGCTACTGAGG
CTGAGGCAGGAGTATTCCTTGAGTCCAGGAGGTCATGGCTGCAGTGAGTTATGATTGTGCCGTCATACTCCAGCCTGAAC
AAAAGAGTGAGACCTTGTCCCTCCCCGCCAAAACCAAACCAAAACAAAACAAAACAAAAAAAAAACACCTAAAAACCCCA
GTGTTTACAGT
"""
    expected = out
    assert observed.strip() == expected.strip()


def test_main_marker_region_mode(capsys):
    arglist = ["marker", "--region", "chr15"]
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split("\n")
    print(out)
    assert len(outlines) == 26 + 1  # markers + 1 header line


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
      Name  NumVars  Extent Chrom     Start       End    Ae   Source
mh06KK-101        2     187  chr6 170280715 170280901 2.196 Kidd2018
mh15KK-058        3     303 chr15  28120285  28120587 3.343 Kidd2018
mh20KK-035        2      31 chr20   2088699   2088729 2.695 Kidd2018
"""
    print(terminal.out)
    assert testout.strip() == terminal.out.strip()


def test_main_marker_panel_plus_query(capsys):
    with NamedTemporaryFile() as panelfile:
        with open(panelfile.name, "w") as fh:
            for marker in ["mh15KK-058", "mh06KK-101", "mh20KK-035"]:
                print(marker, file=fh)
        arglist = ["marker", "--panel", panelfile.name, "--query", "NumVars > 2"]
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


def test_main_marker_bad_code():
    arglist = ["marker", "mh05USC-5pB", "--columns=nausea"]
    args = get_parser().parse_args(arglist)
    with pytest.raises(ValueError, match=r"unsupported format code 'u'"):
        microhapdb.cli.main(args)


@pytest.mark.parametrize(
    "pop,marker,allele,numrows",
    [
        ("--population=Swedish", None, None, 187),
        ("--population=SA000009J", "--marker=mh13KK-218.v1", None, 15),
        (None, "--marker=mh13KK-218.v1", "--allele=C|T|C|T", 103),
        (None, "--marker=mh14PK-72639", None, 86),
        (None, None, None, 124911),
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
    print(out)
    expected = """
         Name  NumVars  Extent Chrom  Start    End                   Positions                 Positions37                                     RSIDs              Source    Ae
mh09KK-033.v1        3      78  chr9 680714 680791        680714;680763;680791        680714;680763;680791           rs10815466;rs9408671;rs17431629 Kidd2018;Turchi2019 3.037
mh09KK-033.v2        4      78  chr9 680714 680791 680714;680763;680768;680791 680714;680763;680768;680791 rs10815466;rs9408671;rs9408672;rs17431629         Staadig2021 3.066
"""
    observed = out
    assert observed.strip() == expected.strip()


def test_ae_pop(capsys):
    arglist = ["marker", "--region=chr18:1-25000000", "--ae-pop=EAS"]
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    exp_out = """
         Name  NumVars  Extent Chrom    Start      End    Ae              Source
  mh18SHY-001       11     181 chr18  1952477  1952657 4.127              Wu2021
 mh18PK-87558        6      47 chr18  1960543  1960589 2.217      vanderGaag2018
 mh18USC-18pA        5      77 chr18  5280018  5280094 1.923      delaPuente2020
   mh18CP-005        4      44 chr18  8892865  8892908 3.620            Kidd2018
  mh18ZBF-002        4      77 chr18 11900703 11900779 3.264             Jin2020
  mh18ZHA-004        4     116 chr18 14315932 14316047 3.017             Sun2020
mh18KK-285.v1        4     136 chr18 24557355 24557490 2.584 Kidd2018;Turchi2019
mh18KK-285.v2        3      59 chr18 24557432 24557490 2.577         Staadig2021
"""
    obs_out = terminal.out
    print(obs_out)
    assert exp_out.strip() == obs_out.strip()


def test_ae_pop_bad_pop():
    arglist = ["marker", "--ae-pop=ABC", "mh18USC-18pA"]
    args = get_parser().parse_args(arglist)
    with pytest.raises(ValueError, match=r'no Ae data for population "ABC"'):
        microhapdb.cli.main(args)
    microhapdb.set_ae_population("1KGP")


def test_marker_offsets_cli(capsys):
    arglist = [
        "marker",
        "--format=offsets",
        "--delta=25",
        "--min-length=200",
        "mh03KK-150.v2",
        "mh11KK-180.v1",
        "mh13KK-217.v1",
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


def test_mhpl8r(capsys):
    arglist = ["frequency", "--marker", "mh02USC-2pA", "--population", "EAS", "--format", "mhpl8r"]
    args = microhapdb.cli.get_parser().parse_args(arglist)
    microhapdb.cli.frequency.main(args)
    terminal = capsys.readouterr()
    result = pandas.read_csv(StringIO(terminal.out), sep="\t")
    print(result)
    assert result.shape == (4, 3)
    assert result.Haplotype.iloc[0] == "A|A|G|A"
    assert result.Frequency.iloc[0] == pytest.approx(0.00780)


def test_mhpl8r_panel(capsys):
    with NamedTemporaryFile(mode="wt") as tempfile:
        print("mh02USC-2pA", "mh08USC-8qA", "mh17USC-17qA", sep="\n", file=tempfile, flush=True)
        arglist = [
            "frequency",
            "--panel",
            tempfile.name,
            "--population",
            "EUR",
            "--format",
            "mhpl8r",
        ]
        args = microhapdb.cli.get_parser().parse_args(arglist)
        microhapdb.cli.frequency.main(args)
    terminal = capsys.readouterr()
    result = pandas.read_csv(StringIO(terminal.out), sep="\t")
    print(result)
    assert result.shape == (14, 3)
    assert result.Haplotype.iloc[8] == "A|G|T"
    assert result.Frequency.iloc[8] == pytest.approx(0.37738)


def test_mhpl8r_multi_pop():
    arglist = [
        "frequency",
        "--marker",
        "mh02USC-2pA",
        "mh08USC-8qA",
        "mh17USC-17qA",
        "--population",
        "SAS",
        "EUR",
        "AFR",
        "--format",
        "mhpl8r",
    ]
    args = microhapdb.cli.get_parser().parse_args(arglist)
    message = r"frequencies for 3 populations recovered, expected only 1"
    with pytest.warns(UserWarning, match=message):
        microhapdb.cli.frequency.main(args)


def test_efm(capsys):
    arglist = [
        "frequency",
        "--population=EUR",
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
    assert result["Allele"].iloc[3] == "C|T|C"
    assert result["mh01USC-1pD"].iloc[3] == pytest.approx(0.093)
    assert pandas.isna(result["mh15USC-15qA"].iloc[3])
    assert result["mh17USC-17pA"].iloc[3] == pytest.approx(0.022)


def test_efm_multi_pop():
    arglist = [
        "frequency",
        "--population",
        "EAS",
        "SAS",
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
    arglist = ["frequency", "--marker", "mh02USC-2pA", "--population", "EAS", "--format", "detail"]
    args = microhapdb.cli.get_parser().parse_args(arglist)
    with pytest.raises(NotImplementedError):
        microhapdb.cli.frequency.main(args)
    args.format = "BoGuS"
    with pytest.raises(ValueError, match=r'unsupported view format "BoGuS"'):
        microhapdb.cli.frequency.main(args)


@pytest.mark.parametrize(
    "arglist",
    [
        ("marker", "BogusMarkerID"),
        ("population", "Atreides"),
        ("frequency", "--marker=mh02USC-2pA", "--population=XYZ"),
        ("frequency", "--marker=NotARealMarkerID", "--population=CEU"),
    ],
)
def test_cli_bad_ids(arglist, capsys):
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    assert terminal.out == ""


def test_cli_locus_name(capsys):
    args = get_parser().parse_args(["marker", "mh01NH-04"])
    microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    observed = terminal.out
    expected = """
        Name  NumVars  Extent Chrom     Start       End    Ae                           Source
mh01NH-04.v2        4     280  chr1 230684605 230684884 4.006 Kidd2018;Turchi2019;Gandotra2020
mh01NH-04.v3        5     280  chr1 230684605 230684884 4.006                      Pakstis2021
mh01NH-04.v1        3      53  chr1 230684832 230684884 3.544          Hiroaki2015;Staadig2021
    """
    print(observed)
    assert observed.strip() == expected.strip()
