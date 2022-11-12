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
    assert Marker.standardize_ids(["SI664625D"]) == ["mh17KK-110"]
    assert Marker.standardize_ids(["MHDBM-19b49670"]) == ["mh11KK-187"]
    assert Marker.standardize_ids(["MHDBM-07c8d144"]) == ["mh04AT-10", "mh04KK-010"]
    assert Marker.standardize_ids(["rs9925859"]) == ["mh16PK-83544"]
    assert Marker.standardize_ids(["rs1533623"]) == ["mh01AT-02", "mh01KK-205", "mh01KKCS-205"]
    assert Marker.standardize_ids(["rs4697751", "SI664918I", "MHDBM-22a6ab26"]) == [
        "mh04CP-007",
        "mh05KK-123",
        "mh11KK-092",
    ]


def test_assumptions():
    num_markers_per_source = [
        198,  # ALFRED
        15,  # 10.1016/j.fsigen.2018.05.008
        40,  # ISFG2019:P597
        26,  # 10.1016/j.legalmed.2015.06.003
        11 - 1,  # 10.1016/j.fsigen.2019.02.018
        10,  # 10.1016/j.fsigen.2018.05.001
        118,  # 10.1016/j.fsigen.2019.102213
        90,  # 10.1016/j.fsigen.2020.102275
        25,  # 10.1016/j.fsigen.2020.102255
        20,  # 10.1098/rsos.191937
        23,  # 10.1002/elps.201900451
        59,  # 10.1007/s00414-020-02483-x
    ]
    assert len(microhapdb.markers) == sum(num_markers_per_source)


def test_markers():
    """
    >>> from microhapdb import Marker
    >>> marker = Marker.from_id("mh18CP-005")
    >>> print(marker)
    mh18CP-005 (chr18:8892864-8892908)
    >>> marker.offsets
    [8892864, 8892893, 8892896, 8892907]
    >>> print(marker.fasta)
    >mh18CP-005 PermID=MHDBM-a85754d3 GRCh38:chr18:8892864-8892908 variants=18,47,50,61 Xref=SI664898P
    GCAGATGTCCTTATACGCAGTGGTGTTAGTTTTAGAAACTGATTCTACGGGTATGCTTGCTCGTGTGTAAAATTATTCAT
    >>> for marker in Marker.from_query("Source == '10.1016/j.fsigen.2018.05.008'"):
    ...   print(marker)
    mh06PK-24844 (chr6:13861392-13861447)
    mh06PK-25713 (chr6:31196949-31197002)
    mh07PK-38311 (chr7:52677450-52677509)
    mh08PK-46625 (chr8:1194352-1194372)
    mh10PK-62104 (chr10:127392565-127392633)
    mh11PK-62906 (chr11:247979-248036)
    mh11PK-63643 (chr11:34415814-34415851)
    mh14PK-72639 (chr14:32203273-32203324)
    mh15PK-75170 (chr15:24802313-24802380)
    mh16PK-83362 (chr16:77999242-77999292)
    mh16PK-83483 (chr16:84516828-84516899)
    mh16PK-83544 (chr16:85934074-85934138)
    mh18PK-87558 (chr18:1960542-1960589)
    mh21PK-MX1s (chr21:41464744-41464824)
    mh22PK-104638 (chr22:44857882-44857955)
    """
    assert microhapdb.markers.shape == (634, 9)
    result = microhapdb.markers[microhapdb.markers.Chrom == "chr19"]
    assert len(result) == 18
    varids = microhapdb.variantmap[microhapdb.variantmap.Marker.isin(result.Name)].Variant.unique()
    assert len(varids) == 74


def test_marker_detail():
    marker = Marker.from_id("mh09AT-17")
    observed = marker.detail
    expected = """
--------------------------------------------------------------[ MicroHapDB ]----
mh09AT-17    a.k.a MHDBM-b34d6dc3

Marker Definition (GRCh38)
    Marker extent
        - chr9:132987175-132987245 (70 bp)
    Target locus
        - chr9:132987165-132987255 (90 bp)
    Constituent variants
        - chromosome offsets: 132987175,132987204,132987244
        - marker offsets: 0,29,69
        - target offsets: 10,39,79
        - cross-references: rs2073578, rs56256724, rs633153
    Observed haplotypes
        - C,A,C
        - C,A,T
        - C,C,C
        - C,C,T
        - T,A,C
        - T,A,T
        - T,C,T


--[ Core Marker Sequence ]--
>mh09AT-17
CTGTCTCTCTGGGCCTCCTCCTCCTAGGAAGGGCGTGCCCTCCTTGCTCCCTCTGGGCTTCCCAGAAACC


--[ Marker Target Sequence with MH alleles (haplotypes) ]--
          *                            *                                       *
CGGCACGTGGCTGTCTCTCTGGGCCTCCTCCTCCTAGGAAGGGCGTGCCCTCCTTGCTCCCTCTGGGCTTCCCAGAAACCGTGGCTATTG
..........C............................A.......................................C..........
..........C............................A.......................................T..........
..........C............................C.......................................C..........
..........C............................C.......................................T..........
..........T............................A.......................................C..........
..........T............................A.......................................T..........
..........T............................C.......................................T..........
--------------------------------------------------------------------------------
"""
    assert observed.strip() == expected.strip()


def test_marker_detail_long():
    marker = Marker.from_id("mh05KK-170", delta=20, minlen=200, extendmode=-1)
    observed = marker.detail
    expected = """
--------------------------------------------------------------[ MicroHapDB ]----
mh05KK-170    a.k.a MHDBM-22ddbb7d, SI664573F

Marker Definition (GRCh38)
    Marker extent
        - chr5:2447909-2448046 (137 bp)
    Target locus
        - chr5:2447866-2448066 (200 bp)
    Constituent variants
        - chromosome offsets: 2447909,2447937,2448031,2448045
        - marker offsets: 0,28,122,136
        - target offsets: 43,71,165,179
        - cross-references: rs370672, rs438055, rs6555108, rs74865590
    Observed haplotypes
        - C,A,A,A
        - C,A,A,G
        - C,A,G,A
        - C,A,G,G
        - C,G,A,A
        - C,G,A,G
        - C,G,G,A
        - C,G,G,G
        - T,A,A,A
        - T,A,A,G
        - T,A,G,A
        - T,A,G,G
        - T,G,A,A
        - T,G,A,G
        - T,G,G,G


--[ Core Marker Sequence ]--
>mh05KK-170
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
mh16PK-83544    a.k.a MHDBM-c80956b4

Marker Definition (GRCh38)
    Marker extent
        - chr16:85934074-85934138 (64 bp)
    Target locus
        - chr16:85934074-85934138 (64 bp)
    Constituent variants
        - chromosome offsets: 85934074,85934082,85934105,85934115,85934128,85934137
        - marker offsets: 0,8,31,41,54,63
        - target offsets: 0,8,31,41,54,63
        - cross-references: rs528179479, rs66509440, rs66804793, rs74032085, rs9925859, rs9929895
    Observed haplotypes
        - C,G,G,C,A,T
        - C,G,G,G,A,T
        - C,G,G,G,G,T
        - T,A,C,G,T,C
        - T,A,G,G,A,T
        - T,A,G,G,T,C


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
mh06PK-24844    a.k.a MHDBM-aa39cbba

Marker Definition (GRCh38)
    Marker extent
        - chr6:13861392-13861447 (55 bp)
    Target locus
        - chr6:13861379-13861460 (81 bp)
    Constituent variants
        - chromosome offsets: 13861392,13861399,13861414,13861421,13861430,13861434,13861438,13861439,13861440,13861446
        - marker offsets: 0,7,22,29,38,42,46,47,48,54
        - target offsets: 13,20,35,42,51,55,59,60,61,67
        - cross-references: rs1196416099, rs1204206, rs1204207, rs1204208, rs34901968, rs35198802, rs376614501, rs545720382, rs546942508, rs553417439
    Observed haplotypes
        - C,C,G,C,C,C,A,A,A,A
        - C,C,G,C,C,C,A,A,G,A
        - G,C,G,T,C,C,A,G,G,A
        - G,C,G,T,G,C,A,G,G,A
        - G,CT,G,C,C,C,A,G,G,A
        - T,C,G,C,C,C,A,A,G,A
        - T,C,G,C,C,T,A,A,G,G
        - T,C,G,C,C,T,T,A,G,G
        - T,C,GA,C,C,T,A,A,G,G


--[ Core Marker Sequence ]--
>mh06PK-24844
TTACATCCAAACGTGAGCAGGAGGAAACTCGGAACATACTGTTTTTAAGAACTAG


--[ Marker Target Sequence with MH alleles (haplotypes) ]--
             *      **              **      *        *   *   ***     *
AGGAAGAAAGTGATTACATCC-AAACGTGAGCAGGAG-GAAACTCGGAACATACTGTTTTTAAGAACTAGTATCACTAGAGTT
.............C......C-..............G-......C........C...C...AAA.....A.............
.............C......C-..............G-......C........C...C...AAG.....A.............
.............G......C-..............G-......T........C...C...AGG.....A.............
.............G......C-..............G-......T........G...C...AGG.....A.............
.............G......CT..............G-......C........C...C...AGG.....A.............
.............T......C-..............G-......C........C...C...AAG.....A.............
.............T......C-..............G-......C........C...T...AAG.....G.............
.............T......C-..............G-......C........C...T...TAG.....G.............
.............T......C-..............GA......C........C...T...AAG.....G.............
--------------------------------------------------------------------------------

--------------------------------------------------------------[ MicroHapDB ]----
mh18PK-87558    a.k.a MHDBM-1e5374f1

Marker Definition (GRCh38)
    Marker extent
        - chr18:1960542-1960589 (47 bp)
    Target locus
        - chr18:1960525-1960606 (81 bp)
    Constituent variants
        - chromosome offsets: 1960542,1960557,1960561,1960566,1960582,1960588
        - marker offsets: 0,15,19,24,40,46
        - target offsets: 17,32,36,41,57,63
        - cross-references: rs28612163, rs28695806, rs62081065, rs62081066, rs9947384, rs9962474
    Observed haplotypes
        - C,A,A,G,A,T
        - C,A,G,C,A,C
        - C,A,G,C,A,T
        - C,A,G,C,T,C
        - C,G,A,G,A,T
        - T,A,G,C,T,C
        - T,G,A,G,A,T


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
mh22PK-104638    a.k.a MHDBM-eb63e78f

Marker Definition (GRCh38)
    Marker extent
        - chr22:44857882-44857955 (73 bp)
    Target locus
        - chr22:44857872-44857965 (93 bp)
    Constituent variants
        - chromosome offsets: 44857882,44857883,44857884,44857891,44857892,44857893,44857907,44857930,44857946,44857949,44857950,44857954
        - marker offsets: 0,1,2,9,10,11,25,48,64,67,68,72
        - target offsets: 10,11,12,19,20,21,35,58,74,77,78,82
        - cross-references: rs10685889, rs10685890, rs113141650, rs117862404, rs62232223, rs62232224, rs62232225, rs62232226, rs71328677, rs7291353
    Observed haplotypes
        - C,A,G,C,G,C,CCTGCC,TTCTT,GTGAG,C,T,G
        - C,A,G,C,G,T,CCTGCC,TTCTT,GTGAG,C,C,T
        - C,G,C,C,G,C,CCTGCC,T,GTGAG,C,T,G
        - C,G,G,C,A,C,CCTGCC,T,G,C,T,G
        - C,G,G,C,G,C,CCTGCC,T,G,C,T,G
        - C,G,G,C,G,C,CCTGCC,T,GTGAG,C,T,G
        - C,G,G,C,G,C,CCTGCC,T,GTGAG,G,T,G
        - C,G,G,C,G,T,C,TTCTT,GTGAG,C,C,T
        - C,G,G,C,G,T,CCTGCC,TTCTT,GTGAG,C,C,T
        - C,G,G,T,G,C,CCTGCC,TTCTT,GTGAG,C,T,G
        - T,G,G,C,G,C,CCTGCC,TTCTT,GTGAG,C,T,G


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
..........CGG......CGC.............CCTGCC.................T----...............G----..CT...G..........
..........CGG......CGC.............CCTGCC.................T----...............GTGAG..CT...G..........
..........CGG......CGC.............CCTGCC.................T----...............GTGAG..GT...G..........
..........CGG......CGT.............C-----.................TTCTT...............GTGAG..CC...T..........
..........CGG......CGT.............CCTGCC.................TTCTT...............GTGAG..CC...T..........
..........CGG......TGC.............CCTGCC.................TTCTT...............GTGAG..CT...G..........
..........TGG......CGC.............CCTGCC.................TTCTT...............GTGAG..CT...G..........
--------------------------------------------------------------------------------
"""
    assert observed.strip() == expected.strip()


def test_marker_fasta():
    marker = Marker.from_id("mh14PK-72639", minlen=70)
    observed = marker.fasta
    expected = """
>mh14PK-72639 PermID=MHDBM-5548fd05 GRCh38:chr14:32203273-32203324 variants=10,11,14,15,28,39,46,57,60
GAGCACGTTGAGAAGAAGTCACCCACAGGCCTTATAAGGCACGGGAGTCTTTCACCTATTCTACTTACGGT
"""
    assert observed.strip() == expected.strip()


def test_marker_fasta_multi(capsys):
    markerids = ["mh05CP-010", "mh13KK-223", "mh14PK-72639"]
    for marker in Marker.from_ids(markerids, delta=15, minlen=150):
        print(marker.fasta)
    observed = capsys.readouterr().out
    expected = """
>mh05CP-010 PermID=MHDBM-df4e5f79 GRCh38:chr5:17903164-17903214 variants=50,69,94,99 Xref=SI664883J
CTGCTCAGAGTTTACATCAGAGTACTTGATGTAAATTACATCAGAGTACGCTGATGTAAATTACATCAGCGTACGCTGAT
GTAAATTACATCAGCGTACTCTGATGTAAATTACATCAGCGTACTCTGATGTAATTTCAGTTTTCTTAAA
>mh13KK-223 PermID=MHDBM-4088f2e4 GRCh38:chr13:110154351-110154505 variants=15,58,75,168 Xref=SI664608E
TTCAGTTGGCTTTTGTGGGAAAGGGAAGCCCTGGGGCTAGGAGAGCAGTCCTTGCCCTCTGGGAAGGGTCCCAGGCGGCA
CTGCCCCAGGAGGGCCTTCGTGGAGGCCACGGCCAGCCCTCGGGTGTTCTCTCCCTAACTCAAGCTTCTGCTTTCAAGCT
CGTGCATGTTGTAGTAGAATGTGT
>mh14PK-72639 PermID=MHDBM-5548fd05 GRCh38:chr14:32203273-32203324 variants=50,51,54,55,68,79,86,97,100
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
    assert marker.slug == "chr11:113425551-113425564"
    assert marker.xrefs == ["MHDBM-f67cb0d4", "SI664596K"]


@pytest.mark.parametrize(
    "markername,slug,length,source",
    [
        ("mh04KK-074", "chr4:55457429-55457526", 97, "ALFRED"),
        ("mh01CP-010", "chr1:85240117-85240140", 23, "10.1016/j.fsigen.2019.02.018"),
        ("mh09KKCS-153", "chr9:101207359-101207606", 247, "10.1016/j.fsigen.2020.102275"),
        ("mh10NH-14", "chr10:11288561-11288613", 52, "10.1016/j.legalmed.2015.06.003"),
        ("mh02ZBF-001", "chr2:99401807-99401896", 89, "10.1002/elps.201900451"),
        ("mh17ZHA-001", "chr17:390129-390249", 120, "10.1098/rsos.191937"),
        ("mh16USC-16pB", "chr16:24314926-24314937", 11, "10.1016/j.fsigen.2019.102213"),
        ("mh13AT-26", "chr13:23191461-23191542", 81, "ISFG2019:P597"),
        ("mh08ZHA-003", "chr8:2914705-2915053", 348, "10.1016/j.fsigen.2020.102255"),
        ("mh15PK-75170", "chr15:24802313-24802380", 67, "10.1016/j.fsigen.2018.05.008"),
        ("mh03LV-06", "chr3:11914400-11914598", 198, "10.1016/j.fsigen.2018.05.001"),
        ("mh15SHY-003", "chr15:92605652-92605846", 194, "10.1007/s00414-020-02483-x"),
    ],
)
def test_all_sources(markername, slug, length, source):
    marker = Marker.from_id(markername)
    assert marker.slug == slug
    assert len(marker) == length
    assert marker.data.Source == source


def test_set_reference():
    coords37 = [
        "22137318,22137395,22137411,22137453",
        "22137395,22137411,22137453",
        "23068395,23068425,23068433",
        "24223721,24223752",
    ]
    coords38 = [
        "24557354,24557431,24557447,24557489",
        "24557431,24557447,24557489",
    ]
    # Make sure they can be swapped in & out multiple times without issues
    for _ in range(4):
        microhapdb.set_reference(37)
        result = Marker.table_from_region("chr18:20000000-25000000")
        assert result.Offsets.tolist() == coords37
        microhapdb.set_reference(38)
        result = Marker.table_from_region("chr18:20000000-25000000")
        assert result.Offsets.tolist() == coords38


@pytest.mark.parametrize(
    "markername,refr,offsets",
    [
        ("mh01KKCS-117", 37, [204633339, 204633396, 204633461, 204633499, 204633525]),
        ("mh06KKCS-104", 38, [165385361, 165385404, 165385450, 165385473, 165385548]),
        ("mh12KKCS-199", 37, [12229785, 12229848, 12229849, 12229885, 12229951]),
        ("mh12KKCS-199", 38, [12076851, 12076914, 12076915, 12076951, 12077017]),
    ],
)
def test_gandotra_offsets(markername, refr, offsets):
    microhapdb.set_reference(refr)
    marker = Marker.from_id(markername)
    microhapdb.set_reference(38)
    observed = marker.offsets
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
    ids = ["mh03AT-09", "mh11KK-180", "mh13KK-217", "mh07USC-7qC"]
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
    expected = sorted([
        "mh12KK-093",
        "mh12KK-045",
        "mh12KK-042",
        "mh12KKCS-046",
        "mh12KK-046",
        "mh12AT-25",
        "mh12CP-003",
        "mh12KKCS-209",
    ])
    assert observed == expected


def test_marker_ids():
    assert Marker.standardize_ids(["BoGUSid"]) == []
    assert Marker.standardize_ids(["mh15KK-058"]) == ["mh15KK-058"]
    assert Marker.standardize_ids(["SI664549I"]) == ["mh01KK-117"]
    assert Marker.standardize_ids(["MHDBM-3d69621c"]) == ["mh11AT-23", "mh11KK-040"]
