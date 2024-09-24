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

from importlib.resources import files
from io import StringIO
import microhapdb
from microhapdb.cli import get_parser
import pandas
import pytest
from tempfile import NamedTemporaryFile


def data_file(path):
    return files("microhapdb") / "tests" / "data" / path


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
    assert len(outlines) == 1 + 125


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
    assert len(outlines) == 1 + 3053


def test_main_marker_notrunc(capsys):
    args = get_parser().parse_args(["marker", "--notrunc"])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    outlines = out.strip().split("\n")
    assert len(outlines) == 1 + 3053


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
                Name  NumVars  Extent Chrom    Start      End     Ae                              Source
   mh19SCUZJ-0001326       17     337 chr19   389061   389397  4.909                             Zhu2023
       mh19WL-007.v2        6     251 chr19   410610   410860  1.126                            Yu2022G3
       mh19WL-007.v1        5     102 chr19   410759   410860  1.000                   Yu2022G1;Yu2022G2
       mh19WL-007.v3        4      89 chr19   410772   410860  1.000                            Yu2022G4
        mh19USC-19pA        3      38 chr19   561779   561816  2.924                      delaPuente2020
         mh19LS-19pA        5     155 chr19   645780   645934  3.724                              Du2023
   mh19SCUZJ-0007484        9     254 chr19   930658   930911  7.592                             Zhu2023
   mh19SCUZJ-0012120        7     232 chr19  1923847  1924078  1.883                             Zhu2023
          mh19WL-024        6     198 chr19  1955525  1955722  4.977                            Yu2022G1
          mh19WL-025        4     170 chr19  3120344  3120513  4.746                            Yu2022G1
   mh19SCUZJ-0018130        7     297 chr19  3125014  3125310  5.460                             Zhu2023
   mh19SCUZJ-0019425        7     335 chr19  3272090  3272424  9.444                             Zhu2023
   mh19SCUZJ-0022447        4     244 chr19  4195774  4196017  5.725                             Zhu2023
          mh19KK-056        2     201 chr19  4852125  4852325  2.608                            Kidd2018
   mh19SCUZJ-0028962       12     346 chr19  5853211  5853556 14.410                             Zhu2023
   mh19SCUZJ-0039302       10     317 chr19  7570487  7570803  9.407                             Zhu2023
      mh19SHY-001.v1        8     185 chr19  7698913  7699097  6.469                              Wu2021
      mh19SHY-001.v2        5     106 chr19  7698992  7699097  6.442           Yu2022G1;Yu2022G2;Zhu2023
      mh19SHY-001.v3        6     280 chr19  7698992  7699271  9.195                            Yu2022G3
      mh19SHY-001.v4        4      67 chr19  7699031  7699097  4.488                            Yu2022G4
   mh19SCUZJ-0046036        7     248 chr19  8787571  8787818  3.649                             Zhu2023
          mh19WL-026        4      84 chr19  9009766  9009849  3.955                            Yu2022G1
   mh19SCUZJ-0056108        6     214 chr19 11146421 11146634  3.835                             Zhu2023
   mh19SCUZJ-0064306        8     170 chr19 13602099 13602268  5.422                             Zhu2023
          mh19CP-007        3      42 chr19 14310740 14310781  3.271                            Kidd2018
          mh19WL-011        3      96 chr19 15931448 15931543  1.803 Yu2022G1;Yu2022G2;Yu2022G3;Yu2022G4
   mh19SCUZJ-0077038        6     233 chr19 15933829 15934061  1.072                             Zhu2023
        mh19USC-19pB        5      66 chr19 16040865 16040930  3.866                      delaPuente2020
   mh19SCUZJ-0081305        5     260 chr19 17174913 17175172  4.770                             Zhu2023
   mh19SCUZJ-0091391        9     269 chr19 20552758 20553026  7.755                             Zhu2023
         mh19ZHA-006        6      63 chr19 20579863 20579925  4.511                             Sun2020
           mh19NH-23        3      95 chr19 22052724 22052818  2.032                         Hiroaki2015
       mh19KK-299.v1        5     154 chr19 22546698 22546851  4.236                 Kidd2018;Turchi2019
       mh19KK-299.v2        7     154 chr19 22546698 22546851  4.253                        Gandotra2020
       mh19KK-299.v3       10     182 chr19 22546698 22546879  4.253                         Pakstis2021
       mh19KK-299.v4        3      63 chr19 22546749 22546811  3.771                         Staadig2021
   mh19SCUZJ-0119796        6     347 chr19 28336605 28336951 10.495                             Zhu2023
         mh19ZHA-007        4     141 chr19 28397316 28397456  4.435                 Kureshi2020;Sun2020
   mh19SCUZJ-0121709        3     120 chr19 28735778 28735897  4.726                             Zhu2023
   mh19SCUZJ-0129055        6     335 chr19 29597152 29597486  6.674                             Zhu2023
   mh19SCUZJ-0137855        5     350 chr19 32851722 32852071  2.958                             Zhu2023
         mh19LS-19qA        3     105 chr19 32967227 32967331  3.048                              Du2023
     mh19USC-19qA.v4        6      78 chr19 33273743 33273820  5.090                           Zhang2023
     mh19USC-19qA.v1        4      46 chr19 33273772 33273817  3.660                      delaPuente2020
     mh19USC-19qA.v3        5      49 chr19 33273772 33273820  5.084                            Yu2022G4
     mh19USC-19qA.v2        6     185 chr19 33273772 33273956  7.147  Yu2022G1;Yu2022G2;Yu2022G3;Zhu2023
   mh19SCUZJ-0146076        8     318 chr19 35194553 35194870  5.130                             Zhu2023
   mh19SCUZJ-0146592       11     307 chr19 35289948 35290254  6.568                             Zhu2023
          mh19WL-022        4     283 chr19 39678442 39678724  3.975                            Yu2022G3
   mh19SCUZJ-0164190        9     338 chr19 40875429 40875766  9.761                             Zhu2023
          mh19WL-013        3      54 chr19 45710986 45711039  2.903                            Yu2022G2
          mh19WL-018        3      92 chr19 46308509 46308600  4.262          Yu2022G1;Yu2022G3;Yu2022G4
          mh19WL-028        5     186 chr19 47379150 47379335  3.727                            Yu2022G1
   mh19SCUZJ-0195134        5     285 chr19 49897779 49898063  7.500                             Zhu2023
          mh19WL-008        5     193 chr19 50811875 50812067  1.993                   Yu2022G1;Yu2022G2
          mh19KK-301        4      64 chr19 50938488 50938551  2.577                 Kidd2018;Turchi2019
          mh19KK-300        7     182 chr19 50947787 50947968  5.407            Gandotra2020;Pakstis2021
          mh19KK-057        3     115 chr19 51654949 51655063  2.498                 Kidd2018;Turchi2019
   mh19SCUZJ-0212784        9     335 chr19 52762327 52762661  6.084                             Zhu2023
   mh19SCUZJ-0215538        5     331 chr19 52912335 52912665  6.622                             Zhu2023
         mh19ZHA-009        5     178 chr19 53129073 53129250  4.286                 Kureshi2020;Sun2020
          mh19WL-020        5     265 chr19 53132309 53132573  7.506                            Yu2022G3
     mh19USC-19qB.v2        7     106 chr19 53714365 53714470  6.254                            Yu2022G1
     mh19USC-19qB.v1        3      27 chr19 53714388 53714414  4.912                      delaPuente2020
     mh19USC-19qB.v3        6      83 chr19 53714388 53714470  6.008                            Yu2022G4
   mh19SCUZJ-0227114        4     241 chr19 54133057 54133297 10.996                             Zhu2023
   mh19SCUZJ-0227616        9     280 chr19 54221872 54222151  1.328                             Zhu2023
   mh19SCUZJ-0238255        7     303 chr19 55526863 55527165  3.786                             Zhu2023
         mh19SHY-002        9     165 chr19 55588421 55588585  3.544                              Wu2021
          mh19WL-027        4     193 chr19 55796289 55796481  4.544                            Yu2022G1
         mh19LS-19qB        3     103 chr19 55911606 55911708  2.435                              Du2023
   mh19SCUZJ-0244558       11     311 chr19 56207216 56207526  5.247                             Zhu2023
          mh19WL-012       11     195 chr19 57013159 57013353  4.385                   Yu2022G1;Yu2022G2
mh19SCUZJ-0249059.v1       17     273 chr19 57245502 57245774 12.735                             Zhu2023
mh19SCUZJ-0249059.v2       10     185 chr19 57245590 57245774  4.085                             Zhu2023
"""
    obs_out = terminal.out
    print(obs_out)
    assert exp_out.strip() == obs_out.strip()


def test_main_marker_fasta_default_delta(capsys):
    args = get_parser().parse_args(["marker", "--format=fasta", "mh01CP-016", "mh06PK-24844"])
    microhapdb.cli.main(args)
    out, err = capsys.readouterr()
    expected = """
>mh01CP-016 GRCh38:chr1:55558994-55559075 mh01CP-016=18,56,62
TGGCACACAACAAGTGCTTATAATGAAAGCATTAGTGAGTAAAAGAGTGATCCCTGGCTTTGAACTCCCTCTAAGTGTAC
C
>mh06PK-24844 GRCh38:chr6:13861379-13861460 mh06PK-24844=13,20,35,42,51,55,59,60,61,67
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
>mh08PK-46625 GRCh38:chr8:1194237-1194487 mh08PK-46625=115,119,127,134
TGCTGGCAAGTTGGAAACACAGGCTCTGCGATTTTGAGAGTGAACCTGCAAGAGACAAGCAGACGTTGGCAGTGCCGCGT
CCCGGCTGGTGGAGGGAGCCCGGATGCCTGGCAGACAGTCAGTGGTCGGTTGGCGGCCGGCCCACATAAGGGCACCATGC
TCACCGTGTCTAGGCAGAGCTGGAGGCTCCTCCTGCCCAGGGCGGCCTCCAGGTGGGGAGGACGGCAGAGCTTCCCTCAG
TCCCACTTTC
>mh13CP-010 GRCh38:chr13:29217935-29218186 mh13CP-010=109,120,141
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
    assert len(outlines) == 101 + 1  # markers + 1 header line


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
mh06KK-101        2     187  chr6 170280715 170280901 2.120 Kidd2018
mh15KK-058        3     303 chr15  28120285  28120587 3.353 Kidd2018
mh20KK-035        2      31 chr20   2088699   2088729 2.703 Kidd2018
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
        ("symmetric", "mh01USC-1pD=70,92,104", "mh22NH-27=47,82,128"),
        ("5", "mh01USC-1pD=120,142,154", "mh22NH-27=73,108,154"),
        ("3", "mh01USC-1pD=20,42,54", "mh22NH-27=20,55,101"),
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
        (None, "--marker=mh13KK-218.v1", "--allele=C|T|C|T", 129),
        (None, "--marker=mh14PK-72639", None, 227),
        (None, None, None, 944737),
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
mh09KK-033.v1        3      78  chr9 680714 680791        680714;680763;680791        680714;680763;680791           rs10815466;rs9408671;rs17431629 Kidd2018;Turchi2019 3.150
mh09KK-033.v2        4      78  chr9 680714 680791 680714;680763;680768;680791 680714;680763;680768;680791 rs10815466;rs9408671;rs9408672;rs17431629         Staadig2021 3.177
"""
    observed = out
    assert observed.strip() == expected.strip()


def test_ae_pop(capsys):
    arglist = ["marker", "--region=chr18:1-25000000", "--ae-pop=EAS"]
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    exp_out = """
                Name  NumVars  Extent Chrom    Start      End     Ae                              Source
           mh18LW-45       12      81 chr18   108028   108108  1.000                           Zhang2023
          mh18WL-002        6     100 chr18   466283   466382  4.056 Yu2022G1;Yu2022G2;Yu2022G3;Yu2022G4
   mh18SCUZJ-0002652        6     343 chr18   626657   626999  1.738                             Zhu2023
   mh18SCUZJ-0003488        7     206 chr18   722819   723024  4.829                             Zhu2023
         mh18LS-18pA        3      97 chr18  1648913  1649009  2.504                              Du2023
         mh18SHY-001       11     181 chr18  1952477  1952657  4.127                              Wu2021
        mh18PK-87558        6      47 chr18  1960543  1960589  2.217                      vanderGaag2018
   mh18SCUZJ-0009510       15     304 chr18  1962681  1962984  1.000                             Zhu2023
   mh18SCUZJ-0011531        4     294 chr18  2326915  2327208  4.654                             Zhu2023
   mh18SCUZJ-0015713        7     269 chr18  3197574  3197842  9.003                             Zhu2023
          mh18WL-028        4     198 chr18  3246291  3246488  4.185                            Yu2022G1
mh18SCUZJ-0020879.v1       11     270 chr18  4477999  4478268 14.775                             Zhu2023
mh18SCUZJ-0020879.v2        6     323 chr18  4478124  4478446  4.430                             Zhu2023
        mh18USC-18pA        5      77 chr18  5280018  5280094  1.923                      delaPuente2020
          mh18WL-026        3     110 chr18  7382444  7382553  3.135                            Yu2022G1
          mh18CP-005        4      44 chr18  8892865  8892908  3.620                            Kidd2018
   mh18SCUZJ-0040982        6     295 chr18  8935646  8935940  5.154                             Zhu2023
   mh18SCUZJ-0042322        6     305 chr18  9555862  9556166  7.195                             Zhu2023
       mh18WL-009.v2        9      97 chr18 10150919 10151015  1.000                            Yu2022G4
       mh18WL-009.v1       10     200 chr18 10150919 10151118  2.000          Yu2022G1;Yu2022G2;Yu2022G3
mh18SCUZJ-0052717.v2        3     110 chr18 11578129 11578238  3.206                              Du2023
mh18SCUZJ-0052717.v1        6     347 chr18 11578129 11578475  7.319                             Zhu2023
         mh18ZBF-002        4      77 chr18 11900703 11900779  3.264                             Jin2020
   mh18SCUZJ-0058695        5     220 chr18 12976468 12976687  6.680                             Zhu2023
mh18SCUZJ-0064707.v1       10     315 chr18 14298002 14298316  6.209                             Zhu2023
mh18SCUZJ-0064707.v2        7     317 chr18 14298131 14298447  2.146                             Zhu2023
         mh18ZHA-004        4     116 chr18 14315932 14316047  3.017                             Sun2020
   mh18SCUZJ-0070441        6     317 chr18 21906065 21906381 10.511                             Zhu2023
       mh18KK-285.v1        4     136 chr18 24557355 24557490  2.584                 Kidd2018;Turchi2019
       mh18KK-285.v2        3      59 chr18 24557432 24557490  2.577                         Staadig2021
   mh18SCUZJ-0076547        7     269 chr18 24716273 24716541  4.878                             Zhu2023
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
        "mh03KK-150.v3",
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


def test_frequency_bad_format():
    arglist = ["frequency", "--marker=mh02USC-2pA", "--population=EAS"]
    args = microhapdb.cli.get_parser().parse_args(arglist)
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
mh01NH-04.v2        4     280  chr1 230684605 230684884 4.029 Kidd2018;Turchi2019;Gandotra2020
mh01NH-04.v3        5     280  chr1 230684605 230684884 4.029                      Pakstis2021
mh01NH-04.v1        3      53  chr1 230684832 230684884 3.566          Hiroaki2015;Staadig2021
    """
    print(observed)
    assert observed.strip() == expected.strip()


def test_cli_summarize(capsys):
    args = get_parser().parse_args(["summarize"])
    microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    observed = terminal.out
    expected = """
[microhaplotypes]
  - 3053 marker definitions
  - 2413 distinct loci
[frequencies]
  - 59753 haplotypes
  - 125 population groups
  - 944736 total microhap frequencies
"""
    print(observed)
    assert observed.strip() == expected.strip()


def test_cli_fasta_locus_multimarker(capsys):
    arglist = [
        "marker",
        "mh14SHY-003.v1",
        "mh14SHY-003.v2",
        "mh14SHY-003.v3",
        "mh14SHY-003.v4",
        "mh02KK-134.v1",
        "mh02KK-134.v2",
        "mh02KK-134.v3",
        "mh02KK-134.v4",
        "--format=fasta",
    ]
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    observed = terminal.out
    expected = """
>mh02KK-134 GRCh38:chr2:160222879-160223013 mh02KK-134.v1=20,44,59,123 mh02KK-134.v2=20,44,59,65,107,123 mh02KK-134.v3=20,44,59,65,72,87,107,123 mh02KK-134.v4=20,44,59
TACCCTTGGCAGGAACCCTCACTACCTAAGGATGGGCAATGGCTTATGAGTGAGAAACACGGAGCCGTGGGAACTCAGAA
TGACATGCTACCTGGAGATTGTGGTAACGCCCTGTTTTTTTGTGGGCATATCTA
>mh14SHY-003 GRCh38:chr14:57983921-57984213 mh14SHY-003.v1=10,14,16,26,102,108,109,161,192,199 mh14SHY-003.v2=102,108,199,262,281 mh14SHY-003.v3=10,14,102,108,199,262,281 mh14SHY-003.v4=102,108,199
GTAGGAGTGATGTACGGGGCACCTACTTGGGGTTCACATGCTGGCCCCTTTATTGAGTTCATTCTGAATCCAGAAGCTTG
GCAGAGTTCAGCCAGATGGCAGGGTGAGCGCCCTGCCTTCCTGGTAGTCTCTTCTTCTGCAAGGGAATAGGAGGCGTTCA
CCCTCCTTTGTTCAAGAGTCTATTTCTAGGGGCCTATCAGCCCAGGGTCCCTTCTCCAGCTTTCTCAGGAGGCCCCACAT
CATCAGGCAATTAGCTCTCTAGTGGGTATAACTGCTACTGCCACAACCACTG
"""
    print(observed)
    assert observed.strip() == expected.strip()


def test_cli_offsets_locus_multimarker(capsys):
    arglist = ["marker", "mh01KK-212", "mh05KK-170", "--delta=20", "--format=offsets"]
    args = get_parser().parse_args(arglist)
    microhapdb.cli.main(args)
    terminal = capsys.readouterr()
    observed = terminal.out
    with open(data_file("multi2.tsv"), "r") as fh:
        expected = fh.read()
    print(observed)
    assert observed.strip() == expected.strip()
