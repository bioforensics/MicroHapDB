# -------------------------------------------------------------------------------------------------
# Copyright (c) 2019, DHS.
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
from microhapdb import Marker
import pandas
import pytest


def test_standardize_ids():
    assert Marker.standardize_ids(["FakeIdentifier"]) == []
    assert Marker.standardize_ids(["mh11KK-089"]) == ["mh11KK-089"]
    assert Marker.standardize_ids(["rs9925859"]) == ["mh16PK-83544"]
    assert Marker.standardize_ids(["rs1533623"]) == [
        "mh01KK-205.v1",
        "mh01KK-205.v2",
        "mh01KK-205.v3",
        "mh01KK-205.v4",
    ]
    assert Marker.standardize_ids(["mh05KK-123", "rs4697751"]) == [
        "mh04CP-007",
        "mh05KK-123",
    ]


def test_assumptions():
    total_markers_per_source = [
        11,  # Chen2019
        22,  # Fan2022
        90,  # Gandotra2020
        26,  # Hiroaki2015
        23,  # Jin2020
        198,  # Kidd2018
        20,  # Kureshi2020
        90,  # Pakstis2021
        45,  # Staadig2021
        25,  # Sun2020
        89,  # Turchi2019
        10,  # Voskoboinik2018
        59,  # Wu2021
        118,  # delaPuente2020
        15,  # vanderGaag2018
    ]
    redundant_markers_per_source = [
        1,  # Chen2019
        2,  # Gandotra2020
        13,  # Pakstis2021
        8,  # Staadig2021
        84,  # Turchi2019
    ]
    expected_markers = sum(total_markers_per_source) - sum(redundant_markers_per_source)
    observed_markers = len(microhapdb.markers)
    assert observed_markers == expected_markers


def test_markers():
    """
    >>> from microhapdb import Marker
    >>> marker = Marker.from_id("mh18CP-005")
    >>> print(marker)
    mh18CP-005 (chr18:8892865-8892908)
    >>> marker.offsets
    [8892864, 8892893, 8892896, 8892907]
    >>> print(marker.fasta)
    >mh18CP-005 GRCh38:chr18:8892846-8892926 variants=18,47,50,61
    GCAGATGTCCTTATACGCAGTGGTGTTAGTTTTAGAAACTGATTCTACGGGTATGCTTGCTCGTGTGTAAAATTATTCAT
    >>> for marker in Marker.from_query("Source == 'vanderGaag2018'"):
    ...   print(marker)
    mh06PK-24844 (chr6:13861393-13861447)
    mh06PK-25713 (chr6:31196950-31197002)
    mh07PK-38311 (chr7:52677451-52677509)
    mh08PK-46625 (chr8:1194353-1194372)
    mh10PK-62104 (chr10:127392566-127392633)
    mh11PK-62906 (chr11:247980-248036)
    mh11PK-63643 (chr11:34415815-34415851)
    mh14PK-72639 (chr14:32203274-32203324)
    mh15PK-75170 (chr15:24802314-24802380)
    mh16PK-83362 (chr16:77999243-77999292)
    mh16PK-83483 (chr16:84516829-84516899)
    mh16PK-83544 (chr16:85934075-85934138)
    mh18PK-87558 (chr18:1960543-1960589)
    mh21PK-MX1s (chr21:41464745-41464824)
    mh22PK-104638 (chr22:44857883-44857955)
    """
    assert microhapdb.markers.shape == (733, 11)
    result = microhapdb.markers[microhapdb.markers.Chrom == "chr19"]
    assert len(result) == 19
    varids = microhapdb.variantmap[microhapdb.variantmap.Marker.isin(result.Name)].Variant.unique()
    assert len(varids) == 77


def test_marker_detail():
    marker = Marker.from_id("mh09KK-157.v3")
    observed = marker.detail
    expected = """
--------------------------------------------------------------[ MicroHapDB ]----
mh09KK-157.v3

Marker Definition
    Marker extent
        - chr9:132987176-132987245 (70 bp)
    Target locus
        - chr9:132987165-132987255 (90 bp)
    Constituent variants
        - chromosome offsets (GRCh37): 135862562, 135862591, 135862631
        - chromosome offsets (GRCh38): 132987175, 132987204, 132987244
        - marker offsets: 0, 29, 69
        - target offsets: 10, 39, 79
        - cross-references: rs56256724, rs2073578, rs633153
    Observed haplotypes
        - C|A|C
        - C|A|T
        - C|C|C
        - C|C|T
        - T|A|T
        - T|C|T


--[ Core Marker Sequence ]--
>mh09KK-157.v3
CTGTCTCTCTGGGCCTCCTCCTCCTAGGAAGGGCGTGCCCTCCTTGCTCCCTCTGGGCTTCCCAGAAACC


--[ Marker Target Sequence with MH alleles (haplotypes) ]--
          *                            *                                       *
CGGCACGTGGCTGTCTCTCTGGGCCTCCTCCTCCTAGGAAGGGCGTGCCCTCCTTGCTCCCTCTGGGCTTCCCAGAAACCGTGGCTATTG
..........C............................A.......................................C..........
..........C............................A.......................................T..........
..........C............................C.......................................C..........
..........C............................C.......................................T..........
..........T............................A.......................................T..........
..........T............................C.......................................T..........
--------------------------------------------------------------------------------
"""
    assert observed.strip() == expected.strip()


def test_marker_detail_long():
    marker = Marker.from_id("mh05KK-170.v1", delta=20, minlen=200, extendmode=-1)
    observed = marker.detail
    expected = """
--------------------------------------------------------------[ MicroHapDB ]----
mh05KK-170.v1

Marker Definition
    Marker extent
        - chr5:2447910-2448046 (137 bp)
    Target locus
        - chr5:2447866-2448066 (200 bp)
    Constituent variants
        - chromosome offsets (GRCh37): 2448023, 2448051, 2448145, 2448159
        - chromosome offsets (GRCh38): 2447909, 2447937, 2448031, 2448045
        - marker offsets: 0, 28, 122, 136
        - target offsets: 43, 71, 165, 179
        - cross-references: rs74865590, rs438055, rs370672, rs6555108
    Observed haplotypes
        - C|A|A|A
        - C|A|A|G
        - C|A|G|A
        - C|A|G|G
        - C|G|A|A
        - C|G|A|G
        - C|G|G|A
        - C|G|G|G
        - T|A|A|A
        - T|A|A|G
        - T|A|G|A
        - T|A|G|G
        - T|G|A|A
        - T|G|A|G
        - T|G|G|G


--[ Core Marker Sequence ]--
>mh05KK-170.v1
CCACAGTTGAAGAGAGAGAGCATGAGACAGCTTGATCGAAATGGTGAAGCTTTGGAGAGATTTTGCGGGGAATGCACTGG
AGAACCATGAAGATGTGGAAACAAACAAACAAACAAAAAACCGCTTTTGCATCTTCA


--[ Marker Target Sequence with MH alleles (haplotypes) ]--
                                           *                           *                                                                                             *             *
TGGAGGACAAAAGTGAACTTGATTTTCTTAACAAAACTGAAGGCCACAGTTGAAGAGAGAGAGCATGAGACAGCTTGATCGAAATGGTGAAGCTTTGGAGAGATTTTGCGGGGAATGCACTGGAGAACCATGAAGATGTGGAAACAAACAAACAAACAAAAAACCGCTTTTGCATCTTCAGACATCTCACTTGTCATCAC
...........................................C...........................A.............................................................................................A.............A....................
...........................................C...........................A.............................................................................................A.............G....................
...........................................C...........................A.............................................................................................G.............A....................
...........................................C...........................A.............................................................................................G.............G....................
...........................................C...........................G.............................................................................................A.............A....................
...........................................C...........................G.............................................................................................A.............G....................
...........................................C...........................G.............................................................................................G.............A....................
...........................................C...........................G.............................................................................................G.............G....................
...........................................T...........................A.............................................................................................A.............A....................
...........................................T...........................A.............................................................................................A.............G....................
...........................................T...........................A.............................................................................................G.............A....................
...........................................T...........................A.............................................................................................G.............G....................
...........................................T...........................G.............................................................................................A.............A....................
...........................................T...........................G.............................................................................................A.............G....................
...........................................T...........................G.............................................................................................G.............G....................
--------------------------------------------------------------------------------
"""
    assert observed.strip() == expected.strip()


def test_marker_detail_short():
    marker = Marker.from_id("mh16PK-83544", delta=0, minlen=0)
    observed = marker.detail
    expected = """
--------------------------------------------------------------[ MicroHapDB ]----
mh16PK-83544

Marker Definition
    Marker extent
        - chr16:85934075-85934138 (64 bp)
    Target locus
        - chr16:85934074-85934138 (64 bp)
    Constituent variants
        - chromosome offsets (GRCh37): 85967680, 85967688, 85967711, 85967721, 85967734, 85967743
        - chromosome offsets (GRCh38): 85934074, 85934082, 85934105, 85934115, 85934128, 85934137
        - marker offsets: 0, 8, 31, 41, 54, 63
        - target offsets: 0, 8, 31, 41, 54, 63
        - cross-references: rs66509440, rs66804793, rs528179479, rs9925859, rs9929895, rs74032085
    Observed haplotypes
        - C|G|G|C|A|T
        - C|G|G|G|A|T
        - C|G|G|G|G|T
        - T|A|C|G|T|C
        - T|A|G|G|A|T
        - T|A|G|G|T|C


--[ Core Marker Sequence ]--
>mh16PK-83544
CGGGCGGTGTGTGGCTGAAGGGCAAGAGAAAGGAGAGACTGCGGCTGGGAACCCATATGGGAAT


--[ Marker Target Sequence with MH alleles (haplotypes) ]--
*       *                      *         *            *        *
CGGGCGGTGTGTGGCTGAAGGGCAAGAGAAAGGAGAGACTGCGGCTGGGAACCCATATGGGAAT
C.......G......................G.........C............A........T
C.......G......................G.........G............A........T
C.......G......................G.........G............G........T
T.......A......................C.........G............T........C
T.......A......................G.........G............A........T
T.......A......................G.........G............T........C
--------------------------------------------------------------------------------
"""
    assert observed.strip() == expected.strip()


def test_marker_detail_multi(capsys):
    for marker in Marker.from_ids(["mh22PK-104638", "mh18PK-87558", "mh06PK-24844"]):
        print(marker.detail)
    observed = capsys.readouterr().out
    expected = """
--------------------------------------------------------------[ MicroHapDB ]----
mh06PK-24844

Marker Definition
    Marker extent
        - chr6:13861393-13861447 (55 bp)
    Target locus
        - chr6:13861379-13861460 (81 bp)
    Constituent variants
        - chromosome offsets (GRCh37): 13861623, 13861630, 13861645, 13861652, 13861661, 13861665, 13861669, 13861670, 13861671, 13861677
        - chromosome offsets (GRCh38): 13861392, 13861399, 13861414, 13861421, 13861430, 13861434, 13861438, 13861439, 13861440, 13861446
        - marker offsets: 0, 7, 22, 29, 38, 42, 46, 47, 48, 54
        - target offsets: 13, 20, 35, 42, 51, 55, 59, 60, 61, 67
        - cross-references: rs1204206, rs376614501, rs1196416099, rs35198802, rs553417439, rs1204207, rs545720382, rs34901968, rs546942508, rs1204208
    Observed haplotypes
        - C|C|G|C|C|C|A|A|A|A
        - C|C|G|C|C|C|A|A|G|A
        - G|CT|G|C|C|C|A|G|G|A
        - G|C|G|T|C|C|A|G|G|A
        - G|C|G|T|G|C|A|G|G|A
        - T|C|GA|C|C|T|A|A|G|G
        - T|C|G|C|C|C|A|A|G|A
        - T|C|G|C|C|T|A|A|G|G
        - T|C|G|C|C|T|T|A|G|G


--[ Core Marker Sequence ]--
>mh06PK-24844
TTACATCCAAACGTGAGCAGGAGGAAACTCGGAACATACTGTTTTTAAGAACTAG


--[ Marker Target Sequence with MH alleles (haplotypes) ]--
             *      **              **      *        *   *   ***     *
AGGAAGAAAGTGATTACATCC-AAACGTGAGCAGGAG-GAAACTCGGAACATACTGTTTTTAAGAACTAGTATCACTAGAGTT
.............C......C-..............G-......C........C...C...AAA.....A.............
.............C......C-..............G-......C........C...C...AAG.....A.............
.............G......CT..............G-......C........C...C...AGG.....A.............
.............G......C-..............G-......T........C...C...AGG.....A.............
.............G......C-..............G-......T........G...C...AGG.....A.............
.............T......C-..............GA......C........C...T...AAG.....G.............
.............T......C-..............G-......C........C...C...AAG.....A.............
.............T......C-..............G-......C........C...T...AAG.....G.............
.............T......C-..............G-......C........C...T...TAG.....G.............
--------------------------------------------------------------------------------

--------------------------------------------------------------[ MicroHapDB ]----
mh18PK-87558

Marker Definition
    Marker extent
        - chr18:1960543-1960589 (47 bp)
    Target locus
        - chr18:1960525-1960606 (81 bp)
    Constituent variants
        - chromosome offsets (GRCh37): 1960543, 1960558, 1960562, 1960567, 1960583, 1960589
        - chromosome offsets (GRCh38): 1960542, 1960557, 1960561, 1960566, 1960582, 1960588
        - marker offsets: 0, 15, 19, 24, 40, 46
        - target offsets: 17, 32, 36, 41, 57, 63
        - cross-references: rs9962474, rs9947384, rs62081065, rs28612163, rs62081066, rs28695806
    Observed haplotypes
        - C|A|A|G|A|T
        - C|A|G|C|A|C
        - C|A|G|C|A|T
        - C|A|G|C|T|C
        - C|G|A|G|A|T
        - T|A|G|C|T|C
        - T|G|A|G|A|T


--[ Core Marker Sequence ]--
>mh18PK-87558
TCAGGTGTTAGCAACGAGGATTTAGAAAAAACAGGTACAAATTATTT


--[ Marker Target Sequence with MH alleles (haplotypes) ]--
                 *              *   *    *               *     *
AGCCTAGCCAAGAGCTGTCAGGTGTTAGCAACGAGGATTTAGAAAAAACAGGTACAAATTATTTCATCACCCAGGTAGTGA
.................C..............A...A....G...............A.....T.................
.................C..............A...G....C...............A.....C.................
.................C..............A...G....C...............A.....T.................
.................C..............A...G....C...............T.....C.................
.................C..............G...A....G...............A.....T.................
.................T..............A...G....C...............T.....C.................
.................T..............G...A....G...............A.....T.................
--------------------------------------------------------------------------------

--------------------------------------------------------------[ MicroHapDB ]----
mh22PK-104638

Marker Definition
    Marker extent
        - chr22:44857883-44857955 (73 bp)
    Target locus
        - chr22:44857872-44857965 (93 bp)
    Constituent variants
        - chromosome offsets (GRCh37): 45253762, 45253763, 45253764, 45253771, 45253772, 45253773, 45253787, 45253810, 45253826, 45253829, 45253830, 45253834
        - chromosome offsets (GRCh38): 44857882, 44857883, 44857884, 44857891, 44857892, 44857893, 44857907, 44857930, 44857946, 44857949, 44857950, 44857954
        - marker offsets: 0, 1, 2, 9, 10, 11, 25, 48, 64, 67, 68, 72
        - target offsets: 10, 11, 12, 19, 20, 21, 35, 58, 74, 77, 78, 82
        - cross-references: rs117862404, rs62232223, rs7291353, rs71328677, rs62232224, rs10685889, rs10685890, rs113141650, rs62232225, rs62232226
    Observed haplotypes
        - C|A|G|C|G|C|CCTGCC|TTCTT|GTGAG|C|T|G
        - C|A|G|C|G|T|CCTGCC|TTCTT|GTGAG|C|C|T
        - C|G|C|C|G|C|CCTGCC|T|GTGAG|C|T|G
        - C|G|G|C|A|C|CCTGCC|T|G|C|T|G
        - C|G|G|C|G|C|CCTGCC|T|GTGAG|C|T|G
        - C|G|G|C|G|C|CCTGCC|T|GTGAG|G|T|G
        - C|G|G|C|G|C|CCTGCC|T|G|C|T|G
        - C|G|G|C|G|T|CCTGCC|TTCTT|GTGAG|C|C|T
        - C|G|G|C|G|T|C|TTCTT|GTGAG|C|C|T
        - C|G|G|T|G|C|CCTGCC|TTCTT|GTGAG|C|T|G
        - T|G|G|C|G|C|CCTGCC|TTCTT|GTGAG|C|T|G


--[ Core Marker Sequence ]--
>mh22PK-104638
CGGTTGTGACGCTCAGCTCACCAGTCCTGCCTACTTGCCAGCAGGTATTCTCAGAGGGACCACAGAGCTGAGG


--[ Marker Target Sequence with MH alleles (haplotypes) ]--
          ***      ***             ******                 *****               *****  **   *
GTGACATTGGCGGTTGTGACGCTCAGCTCACCAGTCCTGCCTACTTGCCAGCAGGTATT----CTCAGAGGGACCACAG----AGCTGAGGGTGACCTGCA
..........CAG......CGC.............CCTGCC.................TTCTT...............GTGAG..CT...G..........
..........CAG......CGT.............CCTGCC.................TTCTT...............GTGAG..CC...T..........
..........CGC......CGC.............CCTGCC.................T----...............GTGAG..CT...G..........
..........CGG......CAC.............CCTGCC.................T----...............G----..CT...G..........
..........CGG......CGC.............CCTGCC.................T----...............GTGAG..CT...G..........
..........CGG......CGC.............CCTGCC.................T----...............GTGAG..GT...G..........
..........CGG......CGC.............CCTGCC.................T----...............G----..CT...G..........
..........CGG......CGT.............CCTGCC.................TTCTT...............GTGAG..CC...T..........
..........CGG......CGT.............C-----.................TTCTT...............GTGAG..CC...T..........
..........CGG......TGC.............CCTGCC.................TTCTT...............GTGAG..CT...G..........
..........TGG......CGC.............CCTGCC.................TTCTT...............GTGAG..CT...G..........
--------------------------------------------------------------------------------
"""
    assert observed.strip() == expected.strip()


def test_marker_fasta():
    marker = Marker.from_id("mh14PK-72639", minlen=70)
    observed = marker.fasta
    expected = """
>mh14PK-72639 GRCh38:chr14:32203263-32203334 variants=10,11,14,15,28,39,46,57,60
GAGCACGTTGAGAAGAAGTCACCCACAGGCCTTATAAGGCACGGGAGTCTTTCACCTATTCTACTTACGGT
"""
    assert observed.strip() == expected.strip()


def test_marker_fasta_multi(capsys):
    markerids = ["mh05CP-010", "mh13KK-223.v1", "mh14PK-72639"]
    for marker in Marker.from_ids(markerids, delta=15, minlen=150):
        print(marker.fasta)
    observed = capsys.readouterr().out
    expected = """
>mh05CP-010 GRCh38:chr5:17903114-17903264 variants=50,69,94,99
CTGCTCAGAGTTTACATCAGAGTACTTGATGTAAATTACATCAGAGTACGCTGATGTAAATTACATCAGCGTACGCTGAT
GTAAATTACATCAGCGTACTCTGATGTAAATTACATCAGCGTACTCTGATGTAATTTCAGTTTTCTTAAA
>mh13KK-223.v1 GRCh38:chr13:110154336-110154520 variants=15,58,75,168
TTCAGTTGGCTTTTGTGGGAAAGGGAAGCCCTGGGGCTAGGAGAGCAGTCCTTGCCCTCTGGGAAGGGTCCCAGGCGGCA
CTGCCCCAGGAGGGCCTTCGTGGAGGCCACGGCCAGCCCTCGGGTGTTCTCTCCCTAACTCAAGCTTCTGCTTTCAAGCT
CGTGCATGTTGTAGTAGAATGTGT
>mh14PK-72639 GRCh38:chr14:32203223-32203374 variants=50,51,54,55,68,79,86,97,100
CTCTGTATCGTTCCAATTTTAGTATATGTGCTGCCGAAGCGAGCACGTTGAGAAGAAGTCACCCACAGGCCTTATAAGGC
ACGGGAGTCTTTCACCTATTCTACTTACGGTGACCGAACCGCGCCCTTTCCTGTCCATCTTGGAGCCTTTG
"""
    assert observed.strip() == expected.strip()


def test_marker_object(capsys):
    marker = Marker.from_id("mh11KK-090", delta=10, minlen=60)
    assert marker.local_to_global(10) == 113425537
    assert marker.global_to_local(113425559) == 32
    assert marker.local_to_global(65) is None
    assert marker.global_to_local(113425599) is None
    assert marker.slug == "chr11:113425552-113425564"


@pytest.mark.parametrize(
    "markername,slug,length,source",
    [
        ("mh04KK-074", "chr4:55457430-55457526", 97, "Kidd2018"),
        ("mh01CP-010", "chr1:85240118-85240140", 23, "Chen2019"),
        ("mh09KK-153.v2", "chr9:101207360-101207606", 247, "Gandotra2020"),
        ("mh10NH-14", "chr10:11288562-11288613", 52, "Hiroaki2015"),
        ("mh02ZBF-001", "chr2:99401808-99401896", 89, "Jin2020"),
        ("mh17ZHA-001", "chr17:390130-390249", 120, "Kureshi2020"),
        ("mh16USC-16pB", "chr16:24314927-24314937", 11, "delaPuente2020"),
        ("mh13KK-213.v3", "chr13:23191462-23191496", 35, "Staadig2021"),
        ("mh08ZHA-003", "chr8:2914706-2915053", 348, "Sun2020"),
        ("mh15PK-75170", "chr15:24802314-24802380", 67, "vanderGaag2018"),
        ("mh03LV-06", "chr3:11914401-11914598", 198, "Voskoboinik2018"),
        ("mh15SHY-003", "chr15:92605653-92605846", 194, "Wu2021"),
        ("mh01FHL-009", "chr1:231954505-231954667", 163, "Fan2022"),
    ],
)
def test_all_sources(markername, slug, length, source):
    marker = Marker.from_id(markername)
    assert marker.slug == slug
    assert len(marker) == length
    assert marker.data.Source == source


@pytest.mark.parametrize(
    "markername,attr,offsets",
    [
        ("mh01KK-117.v2", "offsets37", [204633339, 204633396, 204633461, 204633499, 204633525]),
        ("mh06KK-104", "offsets", [165385361, 165385404, 165385450, 165385473, 165385548]),
        ("mh12KK-199.v1", "offsets37", [12229785, 12229848, 12229849, 12229885, 12229951]),
        ("mh12KK-199.v1", "offsets", [12076851, 12076914, 12076915, 12076951, 12077017]),
    ],
)
def test_gandotra_offsets(markername, attr, offsets):
    marker = Marker.from_id(markername)
    observed = getattr(marker, attr)
    expected = offsets
    assert observed == expected


@pytest.mark.parametrize(
    "marker,mode,offsets",
    [
        ("mh01USC-1pD", 0, "variants=70,92,104"),
        ("mh01USC-1pD", -1, "variants=120,142,154"),
        ("mh01USC-1pD", 1, "variants=20,42,54"),
        ("mh22NH-27", 0, "variants=47,82,128"),
        ("mh22NH-27", -1, "variants=73,108,154"),
        ("mh22NH-27", 1, "variants=20,55,101"),
    ],
)
def test_marker_extendmode(marker, mode, offsets):
    marker = Marker.from_id(marker, delta=20, minlen=175, extendmode=mode)
    assert offsets in marker.defline


def test_marker_extendmode_bad():
    with pytest.raises(ValueError, match=r"invalid literal for int"):
        markers = Marker.from_ids(["mh01USC-1pD", "mh22NH-27"], minlen=175, extendmode="NotAnInt")
        for marker in markers:
            pass
    with pytest.raises(TypeError, match=r"argument must be a string"):
        markers = Marker.from_ids(["mh01USC-1pD", "mh22NH-27"], minlen=175, extendmode=None)
        for marker in markers:
            pass


def test_marker_definition_single():
    marker = Marker.from_id("mh07SHY-002")
    definition = marker.definition
    assert definition.shape == (9, 4)
    assert definition.Offset.tolist() == [10, 36, 45, 67, 68, 90, 116, 121, 122]


def test_marker_definition_multi():
    ids = ["mh03KK-150.v2", "mh11KK-180.v1", "mh13KK-217.v1", "mh07USC-7qC"]
    definition = Marker.definitions_from_ids(ids, delta=25, minlen=200)
    assert definition.shape == (15, 4)
    observed = definition.Offset.tolist()
    expected = [85, 114, 66, 95, 122, 123, 134, 25, 145, 203, 218, 25, 65, 179, 217]
    assert observed == expected
    observed = definition.ChromOffset.tolist()
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


def test_from_region():
    with pytest.raises(ValueError, match='cannot parse region "chr7:123-456-789"'):
        Marker.parse_regionstr("chr7:123-456-789")
    assert len(Marker.table_from_region("chrX")) == 11
    assert len(Marker.table_from_region("chrY")) == 0
    markers = list(Marker.from_region("chr12:100000000-200000000"))
    assert len(markers) == 8
    observed = sorted([marker.name for marker in markers])
    expected = sorted(
        [
            "mh12CP-003",
            "mh12KK-042",
            "mh12KK-045",
            "mh12KK-046.v1",
            "mh12KK-046.v2",
            "mh12KK-046.v3",
            "mh12KK-093",
            "mh12KK-209",
        ]
    )
    assert observed == expected


def test_from_id_no_such_marker():
    with pytest.raises(ValueError, match=r"no such marker 'BoGUSid'"):
        Marker.from_id("BoGUSid")


def test_from_id_no_such_marker():
    with pytest.raises(ValueError, match=r"no such marker 'BoGUSid'"):
        Marker.from_id("BoGUSid")
