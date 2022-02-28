# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


from io import StringIO
import microhapdb
from microhapdb.marker import standardize_ids
import pandas
import pytest


def test_standardize_ids():
    assert list(standardize_ids(['FakeIdentifier']).values) == []
    assert list(standardize_ids(['mh11KK-089']).values) == ['mh11KK-089']
    assert list(standardize_ids(['SI664625D']).values) == ['mh17KK-110']
    assert list(standardize_ids(['MHDBM-19b49670']).values) == ['mh11KK-187']
    assert list(standardize_ids(['MHDBM-07c8d144']).values) == ['mh04KK-010', 'mh04AT-10']
    assert list(standardize_ids(['rs9925859']).values) == ['mh16PK-83544']
    assert list(standardize_ids(['rs1533623']).values) == ['mh01KKCS-205', 'mh01KK-205', 'mh01AT-02']
    assert list(standardize_ids(['rs4697751', 'SI664918I', 'MHDBM-22a6ab26']).values) == ['mh04CP-007', 'mh05KK-123', 'mh11KK-092']


def test_assumptions():
    assert len(microhapdb.markers) == 198 + 15 + 40 + 26 + (11 - 1) + 10 + 118 + 90 + 25


def test_markers():
    """
    >>> import microhapdb
    >>> m = microhapdb.markers
    >>> m[m.Name == 'mh18CP-005']
               Name          PermID Reference  Chrom                          Offsets      Ae      In     Fst  Source
    446  mh18CP-005  MHDBM-a85754d3    GRCh38  chr18  8892864,8892893,8892896,8892907  3.6722  0.0904  0.0059  ALFRED
    >>> m[m.Name == 'mh01KK-117']
              Name          PermID Reference Chrom                                  Offsets      Ae      In     Fst  Source
    31  mh01KK-117  MHDBM-39dc025f    GRCh38  chr1  204664211,204664268,204664371,204664397  4.4565  0.1933  0.0472  ALFRED
    >>> m[m.Name == 'mh11PK-63643']
                 Name          PermID Reference  Chrom                                            Offsets  Ae  In  Fst                        Source
    287  mh11PK-63643  MHDBM-c5ce121f    GRCh38  chr11  34415814,34415816,34415818,34415835,34415836,3... NaN NaN  NaN  10.1016/j.fsigen.2018.05.008
    >>> m[m.Name == 'mh02AT-05']
             Name          PermID Reference Chrom                        Offsets      Ae   In    Fst         Source
    73  mh02AT-05  MHDBM-c3feaba8    GRCh38  chr2  160222899,160222923,160222938  4.5544  0.2  0.152  ISFG2019:P597
    """
    m = microhapdb.markers
    vm = microhapdb.variantmap
    assert m.shape == (532, 9)
    result = m[m.Chrom == 'chr19']
    assert len(result) == 14
    varids = vm[vm.Marker.isin(result.Name)].Variant.unique()
    assert len(varids) == 48


def test_marker_table(capsys):
    marker = microhapdb.markers[microhapdb.markers.Name == 'mh04CP-003']
    microhapdb.marker.print_table(marker)
    testout = '''
      Name         PermID Reference Chrom                 Offsets     Ae     In    Fst Source
mh04CP-003 MHDBM-2be52d8b    GRCh38  chr4 4324722,4324735,4324749 2.9151 0.1668 0.0321 ALFRED
'''
    terminal = capsys.readouterr()
    print(terminal.out)
    assert terminal.out.strip() == testout.strip()


def test_marker_table_multi(capsys):
    markers = microhapdb.markers.query('Name.str.contains("PK")', engine='python').head(n=5)
    microhapdb.marker.print_table(markers)
    testout = '''
         Name          PermID Reference  Chrom                                            Offsets      Ae      In     Fst                        Source
 mh06PK-24844  MHDBM-aa39cbba    GRCh38   chr6  13861392,13861399,13861414,13861421,13861430,1...     NaN     NaN     NaN  10.1016/j.fsigen.2018.05.008
 mh06PK-25713  MHDBM-7d00efdc    GRCh38   chr6  31196949,31196961,31196972,31196985,31196992,3...  2.9544  0.0898  0.0681  10.1016/j.fsigen.2018.05.008
 mh07PK-38311  MHDBM-3ae6dc1b    GRCh38   chr7                52677450,52677456,52677462,52677508  3.2638  0.1000  0.1927  10.1016/j.fsigen.2018.05.008
 mh08PK-46625  MHDBM-840756f3    GRCh38   chr8                    1194352,1194356,1194364,1194371  2.3775  0.1132  0.1395  10.1016/j.fsigen.2018.05.008
 mh10PK-62104  MHDBM-5f9c6cab    GRCh38  chr10  127392565,127392577,127392596,127392610,127392...  2.0068  0.0482  0.1032  10.1016/j.fsigen.2018.05.008
'''
    testoutlong = '''
        Name         PermID Reference Chrom                                                                                   Offsets     Ae     In    Fst                       Source
mh06PK-24844 MHDBM-aa39cbba    GRCh38  chr6 13861392,13861399,13861414,13861421,13861430,13861434,13861438,13861439,13861440,13861446    NaN    NaN    NaN 10.1016/j.fsigen.2018.05.008
mh06PK-25713 MHDBM-7d00efdc    GRCh38  chr6                                     31196949,31196961,31196972,31196985,31196992,31197001 2.9544 0.0898 0.0681 10.1016/j.fsigen.2018.05.008
mh07PK-38311 MHDBM-3ae6dc1b    GRCh38  chr7                                                       52677450,52677456,52677462,52677508 3.2638 0.1000 0.1927 10.1016/j.fsigen.2018.05.008
mh08PK-46625 MHDBM-840756f3    GRCh38  chr8                                                           1194352,1194356,1194364,1194371 2.3775 0.1132 0.1395 10.1016/j.fsigen.2018.05.008
mh10PK-62104 MHDBM-5f9c6cab    GRCh38 chr10                     127392565,127392577,127392596,127392610,127392611,127392620,127392632 2.0068 0.0482 0.1032 10.1016/j.fsigen.2018.05.008
'''
    terminal = capsys.readouterr()
    print(terminal.out)
    assert terminal.out.strip() == testout.strip() or terminal.out.strip() == testoutlong.strip()


def test_marker_table_multi_notrunc(capsys):
    markers = microhapdb.markers.query('Name.str.contains("PK")', engine='python').head(n=5)
    microhapdb.marker.print_table(markers, trunc=False)
    testoutlong = '''
        Name         PermID Reference Chrom                                                                                   Offsets     Ae     In    Fst                       Source
mh06PK-24844 MHDBM-aa39cbba    GRCh38  chr6 13861392,13861399,13861414,13861421,13861430,13861434,13861438,13861439,13861440,13861446    NaN    NaN    NaN 10.1016/j.fsigen.2018.05.008
mh06PK-25713 MHDBM-7d00efdc    GRCh38  chr6                                     31196949,31196961,31196972,31196985,31196992,31197001 2.9544 0.0898 0.0681 10.1016/j.fsigen.2018.05.008
mh07PK-38311 MHDBM-3ae6dc1b    GRCh38  chr7                                                       52677450,52677456,52677462,52677508 3.2638 0.1000 0.1927 10.1016/j.fsigen.2018.05.008
mh08PK-46625 MHDBM-840756f3    GRCh38  chr8                                                           1194352,1194356,1194364,1194371 2.3775 0.1132 0.1395 10.1016/j.fsigen.2018.05.008
mh10PK-62104 MHDBM-5f9c6cab    GRCh38 chr10                     127392565,127392577,127392596,127392610,127392611,127392620,127392632 2.0068 0.0482 0.1032 10.1016/j.fsigen.2018.05.008
'''
    terminal = capsys.readouterr()
    print(terminal.out)
    assert terminal.out.strip() == testoutlong.strip()


def test_marker_detail(capsys):
    marker = microhapdb.markers[microhapdb.markers.Name == 'mh09AT-17']
    microhapdb.marker.print_detail(marker, delta=10, minlen=90)
    testout = '''
--------------------------------------------------------------[ MicroHapDB ]----
mh09AT-17    a.k.a MHDBM-b34d6dc3

Marker Definition (GRCh38)
    Marker extent
        - chr9:132987175-132987245 (70 bp)
    Target amplicon
        - chr9:132987165-132987255 (90 bp)
    Constituent variants
        - chromosome offsets: 132987175,132987204,132987244
        - marker offsets: 0,29,69
        - amplicon offsets: 10,39,79
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


--[ Target Amplicon Sequence with Haplotypes ]--
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
'''
    terminal = capsys.readouterr()
    assert terminal.out.strip() == testout.strip()


def test_marker_detail_long(capsys):
    marker = microhapdb.markers[microhapdb.markers.Name == 'mh05KK-170']
    microhapdb.marker.print_detail(marker, delta=20, minlen=125)
    testout = '''
--------------------------------------------------------------[ MicroHapDB ]----
mh05KK-170    a.k.a MHDBM-22ddbb7d, SI664573F

Marker Definition (GRCh38)
    Marker extent
        - chr5:2447909-2448046 (137 bp)
    Target amplicon
        - chr5:2447889-2448066 (177 bp)
    Constituent variants
        - chromosome offsets: 2447909,2447937,2448031,2448045
        - marker offsets: 0,28,122,136
        - amplicon offsets: 20,48,142,156
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


--[ Target Amplicon Sequence with Haplotypes ]--
                    *                           *                                                                                             *             *
TTTCTTAACAAAACTGAAGGCCACAGTTGAAGAGAGAGAGCATGAGACAGCTTGATCGAAATGGTGAAGCTTTGGAGAGATTTTGCGGGGAATGCACTGGAGAACCATGAAGATGTGGAAACAAACAAACAAACAAAAAACCGCTTTTGCATCTTCAGACATCTCACTTGTCATCAC
....................C...........................A.............................................................................................A.............A....................
....................C...........................A.............................................................................................A.............G....................
....................C...........................A.............................................................................................G.............A....................
....................C...........................A.............................................................................................G.............G....................
....................C...........................G.............................................................................................A.............A....................
....................C...........................G.............................................................................................A.............G....................
....................C...........................G.............................................................................................G.............A....................
....................C...........................G.............................................................................................G.............G....................
....................T...........................A.............................................................................................A.............A....................
....................T...........................A.............................................................................................A.............G....................
....................T...........................A.............................................................................................G.............A....................
....................T...........................A.............................................................................................G.............G....................
....................T...........................G.............................................................................................A.............A....................
....................T...........................G.............................................................................................A.............G....................
....................T...........................G.............................................................................................G.............G....................
--------------------------------------------------------------------------------
'''
    terminal = capsys.readouterr()
    assert terminal.out.strip() == testout.strip()


def test_marker_detail_short(capsys):
    marker = microhapdb.markers[microhapdb.markers.Name == 'mh16PK-83544']
    microhapdb.marker.print_detail(marker, delta=0, minlen=0)
    testout = '''
--------------------------------------------------------------[ MicroHapDB ]----
mh16PK-83544    a.k.a MHDBM-c80956b4

Marker Definition (GRCh38)
    Marker extent
        - chr16:85934074-85934138 (64 bp)
    Target amplicon
        - chr16:85934074-85934138 (64 bp)
    Constituent variants
        - chromosome offsets: 85934074,85934082,85934105,85934115,85934128,85934137
        - marker offsets: 0,8,31,41,54,63
        - amplicon offsets: 0,8,31,41,54,63
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


--[ Target Amplicon Sequence with Haplotypes ]--
*       *                      *         *            *        *
CGGGCGGTGTGTGGCTGAAGGGCAAGAGAAAGGAGAGACTGCGGCTGGGAACCCATATGGGAAT
C.......G......................G.........C............A........T
C.......G......................G.........G............A........T
C.......G......................G.........G............G........T
T.......A......................C.........G............T........C
T.......A......................G.........G............A........T
T.......A......................G.........G............T........C
--------------------------------------------------------------------------------


'''
    terminal = capsys.readouterr()
    assert terminal.out.strip() == testout.strip()


def test_marker_detail_multi(capsys):
    markerids = ['mh22PK-104638', 'mh18PK-87558', 'mh06PK-24844']
    markers = microhapdb.markers[microhapdb.markers.Name.isin(markerids)]
    microhapdb.marker.print_detail(markers, delta=10, minlen=80)
    testout = '''
--------------------------------------------------------------[ MicroHapDB ]----
mh06PK-24844    a.k.a MHDBM-aa39cbba

Marker Definition (GRCh38)
    Marker extent
        - chr6:13861392-13861447 (55 bp)
    Target amplicon
        - chr6:13861379-13861460 (81 bp)
    Constituent variants
        - chromosome offsets: 13861392,13861399,13861414,13861421,13861430,13861434,13861438,13861439,13861440,13861446
        - marker offsets: 0,7,22,29,38,42,46,47,48,54
        - amplicon offsets: 13,20,35,42,51,55,59,60,61,67
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


--[ Target Amplicon Sequence with Haplotypes ]--
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
    Target amplicon
        - chr18:1960525-1960606 (81 bp)
    Constituent variants
        - chromosome offsets: 1960542,1960557,1960561,1960566,1960582,1960588
        - marker offsets: 0,15,19,24,40,46
        - amplicon offsets: 17,32,36,41,57,63
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


--[ Target Amplicon Sequence with Haplotypes ]--
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
    Target amplicon
        - chr22:44857872-44857965 (93 bp)
    Constituent variants
        - chromosome offsets: 44857882,44857883,44857884,44857891,44857892,44857893,44857907,44857930,44857946,44857949,44857950,44857954
        - marker offsets: 0,1,2,9,10,11,25,48,64,67,68,72
        - amplicon offsets: 10,11,12,19,20,21,35,58,74,77,78,82
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


--[ Target Amplicon Sequence with Haplotypes ]--
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
'''
    terminal = capsys.readouterr()
    assert terminal.out.strip() == testout.strip()


def test_marker_fasta(capsys):
    markers = microhapdb.markers[microhapdb.markers.Name == 'mh14PK-72639']
    microhapdb.marker.print_fasta(markers, delta=10, minlen=70)
    testout = '''
>mh14PK-72639 PermID=MHDBM-5548fd05 GRCh38:chr14:32203273-32203323 variants=10,11,14,15,28,39,46,57,60
GAGCACGTTGAGAAGAAGTCACCCACAGGCCTTATAAGGCACGGGAGTCTTTCACCTATTCTACTTACGGT
'''
    terminal = capsys.readouterr()
    print(terminal.out)
    assert terminal.out.strip() == testout.strip()


def test_marker_fasta_multi(capsys):
    markerids = ['mh05CP-010', 'mh13KK-223', 'mh14PK-72639']
    markers = microhapdb.markers[microhapdb.markers.Name.isin(markerids)]
    microhapdb.marker.print_fasta(markers, delta=15, minlen=150)
    testout = '''
>mh05CP-010 PermID=MHDBM-df4e5f79 GRCh38:chr5:17903164-17903213 variants=50,69,94,99 Xref=SI664883J
CTGCTCAGAGTTTACATCAGAGTACTTGATGTAAATTACATCAGAGTACGCTGATGTAAATTACATCAGCGTACGCTGAT
GTAAATTACATCAGCGTACTCTGATGTAAATTACATCAGCGTACTCTGATGTAATTTCAGTTTTCTTAAA
>mh13KK-223 PermID=MHDBM-4088f2e4 GRCh38:chr13:110154351-110154504 variants=15,58,75,168 Xref=SI664608E
TTCAGTTGGCTTTTGTGGGAAAGGGAAGCCCTGGGGCTAGGAGAGCAGTCCTTGCCCTCTGGGAAGGGTCCCAGGCGGCA
CTGCCCCAGGAGGGCCTTCGTGGAGGCCACGGCCAGCCCTCGGGTGTTCTCTCCCTAACTCAAGCTTCTGCTTTCAAGCT
CGTGCATGTTGTAGTAGAATGTGT
>mh14PK-72639 PermID=MHDBM-5548fd05 GRCh38:chr14:32203273-32203323 variants=50,51,54,55,68,79,86,97,100
CTCTGTATCGTTCCAATTTTAGTATATGTGCTGCCGAAGCGAGCACGTTGAGAAGAAGTCACCCACAGGCCTTATAAGGC
ACGGGAGTCTTTCACCTATTCTACTTACGGTGACCGAACCGCGCCCTTTCCTGTCCATCTTGGAGCCTTTG
'''
    terminal = capsys.readouterr()
    print(terminal.out)
    assert terminal.out.strip() == testout.strip()


def test_amplicon_object(capsys):
    amp = microhapdb.marker.TargetAmplicon('mh11KK-090', delta=10, minlen=60)
    assert amp.local_to_global(10) == 113425537
    assert amp.global_to_local(113425559) == 32
    assert amp.local_to_global(65) is None
    assert amp.global_to_local(113425599) is None
    print(amp)
    testout = '''
--------------------------------------------------------------[ MicroHapDB ]----
mh11KK-090    a.k.a MHDBM-f67cb0d4, SI664596K

Marker Definition (GRCh38)
    Marker extent
        - chr11:113425551-113425564 (13 bp)
    Target amplicon
        - chr11:113425527-113425588 (61 bp)
    Constituent variants
        - chromosome offsets: 113425551,113425563
        - marker offsets: 0,12
        - amplicon offsets: 24,36
        - cross-references: rs1079597, rs1079598
    Observed haplotypes
        - A,C
        - A,T
        - G,C
        - G,T


--[ Core Marker Sequence ]--
>mh11KK-090
AGATTCGCCTTTC


--[ Target Amplicon Sequence with Haplotypes ]--
                        *           *
AAGGGCAGCAGGAACCACATGATCAGATTCGCCTTTCGAATAGGTGATTCTGACAGCACTG
........................A...........C........................
........................A...........T........................
........................G...........C........................
........................G...........T........................
--------------------------------------------------------------------------------
'''
    terminal = capsys.readouterr()
    assert terminal.out.strip() == testout.strip()


@pytest.mark.parametrize('name,data', [
    ('mh04KK-010', 'mh04KK-010 MHDBM-07c8d144    GRCh38  chr4 1985210,1985244 2.9226 0.1306 0.0406 ALFRED'),
    ('mh08PK-46625', 'mh08PK-46625 MHDBM-840756f3    GRCh38  chr8 1194352,1194356,1194364,1194371 2.3775 0.1132 0.1395 10.1016/j.fsigen.2018.05.008'),
    ('mh04AT-10', 'mh04AT-10 MHDBM-07c8d144    GRCh38  chr4 1985210,1985244 2.9226 0.1306 0.0406 ISFG2019:P597'),
    ('mh01NH-03', 'mh01NH-03 MHDBM-e7a95c5e    GRCh38  chr1 184807944,184807966,184808042 1.9085 0.1987 0.1725 10.1016/j.legalmed.2015.06.003'),
    ('mh04CP-004', 'mh04CP-004 MHDBM-8408d717    GRCh38  chr4 7402842,7402854,7402870 2.6744 0.0477 0.061 10.1016/j.fsigen.2019.02.018'),
    ('mh03LV-07', 'mh03LV-07 MHDBM-5f7e29b6    GRCh38  chr3 5783508,5783509,5783518,5783523,5783525,5783531,5783541,5783542,5783543,5783544,5783552,5783562,5783564,5783571,5783577,5783607,5783608,5783611,5783612,5783617,5783618,5783619,5783623,5783626,5783635,5783648,5783652,5783653,5783663,5783664,5783671,5783672,5783673,5783676,5783677,5783678,5783681,5783684,5783687,5783695,5783704,5783705 14.0275 1.081 0.057 10.1016/j.fsigen.2018.05.001'),
    ('mh09USC-9pB', 'mh09USC-9pB MHDBM-7da7af40    GRCh38  chr9 31196676,31196714,31196731,31196744 3.0919 0.1616 0.0574 10.1016/j.fsigen.2019.102213'),
    ('mh13KKCS-223', 'mh13KKCS-223 MHDBM-3ca7e2fc    GRCh38 chr13 110154341,110154351,110154394,110154411,110154438,110154441,110154485,110154504 NaN NaN  NaN 10.1016/j.fsigen.2020.102275')
])
def test_all_sources(name, data, capsys):
    marker = microhapdb.markers[microhapdb.markers.Name == name]
    microhapdb.marker.print_table(marker, trunc=False)
    terminal = capsys.readouterr()
    print(terminal.out)
    assert data in terminal.out


def test_set_reference():
    coords37 = [
        '22137318,22137395,22137411,22137453',
        '22137395,22137411,22137453',
        '23068395,23068425,23068433',
    ]
    coords38 = [
        '24557354,24557431,24557447,24557489',
        '24557431,24557447,24557489',
    ]
    # Make sure they can be swapped in & out multiple times without issues
    for _ in range(4):
        microhapdb.set_reference(37)
        result = microhapdb.retrieve.by_region('chr18:20000000-25000000')
        assert result.Offsets.tolist() == coords37
        microhapdb.set_reference(38)
        result = microhapdb.retrieve.by_region('chr18:20000000-25000000')
        assert result.Offsets.tolist() == coords38


@pytest.mark.parametrize('markername,refr,offsets', [
    (
        'mh11KKCS-180', 37,
        '1690724,1690769,1690790,1690824,1690886,1690910,1690949,1690961,1690968,1690983'
    ),
    (
        'mh11KKCS-180', 38,
        '1669494,1669539,1669560,1669594,1669656,1669680,1669719,1669731,1669738,1669753'
    ),
    (
        'mh12KKCS-199', 37,
        '12229785,12229848,12229849,12229885,12229951'
    ),
    (
        'mh12KKCS-199', 38,
        '12076851,12076914,12076915,12076951,12077017'
    ),
])
def test_gandotra_offsets(markername, refr, offsets, capsys):
    microhapdb.set_reference(refr)
    marker = microhapdb.markers[microhapdb.markers.Name == markername]
    microhapdb.marker.print_detail(marker)
    terminal = capsys.readouterr()
    microhapdb.set_reference(38)
    assert offsets in terminal.out


@pytest.mark.parametrize('markername', ['mh0XUSC-XqD'])
def test_marker_no_freq(markername, capsys):
    marker = microhapdb.markers[microhapdb.markers.Name == markername]
    microhapdb.marker.print_detail(marker)
    terminal = capsys.readouterr()
    message = 'Unable to display a full detail view for markers without frequency information'
    assert message in terminal.err


@pytest.mark.parametrize('mode,offsets1,offsets2', [
    (0, 'variants=70,92,104', 'variants=47,82,128'),
    (-1, 'variants=120,142,154', 'variants=73,108,154'),
    (1, 'variants=20,42,54', 'variants=20,55,101'),
])
def test_marker_extendmode(mode, offsets1, offsets2, capsys):
    markers = microhapdb.markers[microhapdb.markers.Name.isin(['mh01USC-1pD', 'mh22NH-27'])]
    microhapdb.marker.print_fasta(markers, delta=20, minlen=175, extendmode=mode)
    terminal = capsys.readouterr()
    assert offsets1 in terminal.out
    assert offsets2 in terminal.out


def test_marker_extendmode_bad():
    markers = microhapdb.markers[microhapdb.markers.Name.isin(['mh01USC-1pD', 'mh22NH-27'])]
    with pytest.raises(ValueError, match=r'invalid literal for int'):
        microhapdb.marker.print_fasta(markers, delta=20, minlen=175, extendmode="NotAnInt")
    with pytest.raises(TypeError, match=r'argument must be a string'):
        microhapdb.marker.print_fasta(markers, delta=20, minlen=175, extendmode=None)


def test_marker_offsets(capsys):
    ids = ['mh03AT-09', 'mh11KK-180', 'mh13KK-217', 'mh07USC-7qC']
    markers = microhapdb.markers[microhapdb.markers.Name.isin(ids)]
    microhapdb.marker.print_offsets(markers, delta=25, minlen=200)
    terminal = capsys.readouterr()
    result = pandas.read_csv(StringIO(terminal.out), sep='\t')
    assert result.shape == (15, 2)
    observed = list(result.Offset)
    expected = [85, 114, 66, 95, 122, 123, 134, 25, 145, 203, 218, 25, 65, 179, 217]
    assert observed == expected
