# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


import microhapdb
from microhapdb.marker import standardize_ids
import pytest


def test_standardize_ids():
    assert list(standardize_ids(['FakeIdentifier']).values) == []
    assert list(standardize_ids(['mh11KK-089']).values) == ['mh11KK-089']
    assert list(standardize_ids(['SI664625D']).values) == ['mh17KK-110']
    assert list(standardize_ids(['MHDBM-19b49670']).values) == ['mh11KK-187']
    assert list(standardize_ids(['MHDBM-07c8d144']).values) == ['mh04KK-010', 'mh04AT-10']
    assert list(standardize_ids(['rs9925859']).values) == ['mh16PK-83544']
    assert list(standardize_ids(['rs1533623']).values) == ['mh01KK-205', 'mh01AT-02']
    assert list(standardize_ids(['rs4697751', 'SI664918I', 'MHDBM-22a6ab26']).values) == ['mh04CP-007', 'mh05KK-123', 'mh11KK-092']


def test_assumptions():
    assert len(microhapdb.markers) == 198 + 15 + 40 + 26 + 11


def test_markers():
    """
    >>> import microhapdb
    >>> m = microhapdb.markers
    >>> m[m.Name == 'mh18CP-005']
               Name          PermID Reference  Chrom                          Offsets   AvgAe  Source
    250  mh18CP-005  MHDBM-a85754d3    GRCh38  chr18  8892864,8892893,8892896,8892907  3.6635  ALFRED
    >>> m[m.Name == 'mh01KK-117']
              Name          PermID Reference Chrom                                  Offsets  AvgAe  Source
    18  mh01KK-117  MHDBM-39dc025f    GRCh38  chr1  204664211,204664268,204664371,204664397  3.933  ALFRED
    >>> m[m.Name == 'mh11PK-63643']
                 Name          PermID Reference  Chrom                                            Offsets  AvgAe                        Source
    157  mh11PK-63643  MHDBM-c5ce121f    GRCh38  chr11  34415814,34415816,34415818,34415835,34415836,3...  4.033  10.1016/j.fsigen.2018.05.008
    >>> m[m.Name == 'mh02AT-05']
             Name          PermID Reference Chrom                        Offsets   AvgAe         Source
    40  mh02AT-05  MHDBM-c3feaba8    GRCh38  chr2  160222899,160222923,160222938  5.1944  ISFG2019:P597
    """
    m = microhapdb.markers
    vm = microhapdb.variantmap
    assert m.shape == (290, 7)
    result = m[m.Chrom == 'chr19']
    assert len(result) == 7
    varids = vm[vm.Marker.isin(result.Name)].Variant.unique()
    assert len(varids) == 20


def test_marker_table(capsys):
    marker = microhapdb.markers[microhapdb.markers.Name == 'mh04CP-003']
    microhapdb.marker.print_table(marker)
    testout = '''
       Name          PermID Reference Chrom                  Offsets   AvgAe  Source
 mh04CP-003  MHDBM-2be52d8b    GRCh38  chr4  4324722,4324735,4324749  2.9594  ALFRED
'''
    terminal = capsys.readouterr()
    assert terminal.out.strip() == testout.strip()


def test_marker_table_multi(capsys):
    markers = microhapdb.markers.query('Name.str.contains("PK")').head(n=5)
    microhapdb.marker.print_table(markers)
    testout = '''
         Name          PermID Reference  Chrom                                            Offsets   AvgAe                        Source
 mh06PK-24844  MHDBM-aa39cbba    GRCh38   chr6  13861392,13861399,13861414,13861421,13861430,1...  2.4481  10.1016/j.fsigen.2018.05.008
 mh06PK-25713  MHDBM-7d00efdc    GRCh38   chr6  31196949,31196961,31196972,31196985,31196992,3...  3.0005  10.1016/j.fsigen.2018.05.008
 mh07PK-38311  MHDBM-3ae6dc1b    GRCh38   chr7                52677450,52677456,52677462,52677508  2.3053  10.1016/j.fsigen.2018.05.008
 mh08PK-46625  MHDBM-840756f3    GRCh38   chr8                    1194352,1194356,1194364,1194371  2.6273  10.1016/j.fsigen.2018.05.008
 mh10PK-62104  MHDBM-5f9c6cab    GRCh38  chr10  127392565,127392577,127392596,127392610,127392...  2.0199  10.1016/j.fsigen.2018.05.008
'''
    terminal = capsys.readouterr()
    print(terminal.out)
    assert terminal.out.strip() == testout.strip()


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
        - cross-references: rs56256724, rs2073578, rs633153
    Observed alleles
        - C,A,C
        - C,A,T
        - C,C,T
        - T,A,T


--[ Marker Sequence ]--
>mh09AT-17
CTGTCTCTCTGGGCCTCCTCCTCCTAGGAAGGGCGTGCCCTCCTTGCTCCCTCTGGGCTTCCCAGAAACC


--[ Target Amplicon Sequence with Alleles ]--
          *                            *                                       *
CGGCACGTGGCTGTCTCTCTGGGCCTCCTCCTCCTAGGAAGGGCGTGCCCTCCTTGCTCCCTCTGGGCTTCCCAGAAACCGTGGCTATTG
..........C............................A.......................................C..........
..........C............................A.......................................T..........
..........C............................C.......................................T..........
..........T............................A.......................................T..........
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
        - cross-references: rs74865590, rs438055, rs370672, rs6555108
    Observed alleles
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


--[ Marker Sequence ]--
>mh05KK-170
CCACAGTTGAAGAGAGAGAGCATGAGACAGCTTGATCGAAATGGTGAAGCTTTGGAGAGATTTTGCGGGGAATGCACTGG
AGAACCATGAAGATGTGGAAACAAACAAACAAACAAAAAACCGCTTTTGCATCTTCA


--[ Target Amplicon Sequence with Alleles ]--
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
        - cross-references: rs66509440, rs66804793, rs528179479, rs9925859, rs9929895, rs74032085
    Observed alleles
        - C,G,G,C,A,T
        - C,G,G,G,A,T
        - C,G,G,G,G,T
        - T,A,C,G,T,C
        - T,A,G,G,T,C


--[ Marker Sequence ]--
>mh16PK-83544
CGGGCGGTGTGTGGCTGAAGGGCAAGAGAAAGGAGAGACTGCGGCTGGGAACCCATATGGGAAT


--[ Target Amplicon Sequence with Alleles ]--
*       *                      *         *            *        *
CGGGCGGTGTGTGGCTGAAGGGCAAGAGAAAGGAGAGACTGCGGCTGGGAACCCATATGGGAAT
C.......G......................G.........C............A........T
C.......G......................G.........G............A........T
C.......G......................G.........G............G........T
T.......A......................C.........G............T........C
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
        - cross-references: rs1204206, rs1287207523, rs1196416099, rs35198802, rs553417439, rs1204207, rs545720382, rs1427107855, rs546942508, rs1204208
    Observed alleles
        - C,C,G,C,C,C,A,A,A,A
        - C,C,G,C,C,C,A,A,G,A
        - G,C,G,T,C,C,A,G,G,A
        - G,C,G,T,G,C,A,G,G,A
        - G,CT,G,C,C,C,A,G,G,A
        - T,C,G,C,C,C,A,A,G,A
        - T,C,G,C,C,T,A,A,G,G
        - T,C,G,C,C,T,T,A,G,G
        - T,C,GA,C,C,T,A,A,G,G


--[ Marker Sequence ]--
>mh06PK-24844
TTACATCCAAACGTGAGCAGGAGGAAACTCGGAACATACTGTTTTTAAGAACTAG


--[ Target Amplicon Sequence with Alleles ]--
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
        - cross-references: rs9962474, rs9947384, rs62081065, rs28612163, rs62081066, rs28695806
    Observed alleles
        - C,A,A,G,A,T
        - C,A,G,C,A,C
        - C,A,G,C,T,C
        - T,G,A,G,A,T


--[ Marker Sequence ]--
>mh18PK-87558
TCAGGTGTTAGCAACGAGGATTTAGAAAAAACAGGTACAAATTATTT


--[ Target Amplicon Sequence with Alleles ]--
                 *              *   *    *               *     *
AGCCTAGCCAAGAGCTGTCAGGTGTTAGCAACGAGGATTTAGAAAAAACAGGTACAAATTATTTCATCACCCAGGTAGTGA
.................C..............A...A....G...............A.....T.................
.................C..............A...G....C...............A.....C.................
.................C..............A...G....C...............T.....C.................
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
        - cross-references: rs117862404, rs62232223, rs7291353, rs71328677, rs62232224, rs938165705, rs113141650, rs62232225, rs62232226
    Observed alleles
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


--[ Marker Sequence ]--
>mh22PK-104638
CGGTTGTGACGCTCAGCTCACCAGTCCTGCCTACTTGCCAGCAGGTATTCTCAGAGGGACCACAGAGCTGAGG


--[ Target Amplicon Sequence with Alleles ]--
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


def test_amplicon_object(capsys):
    amp = microhapdb.marker.TargetAmplicon('mh11KK-090', delta=10, minlen=60)
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
        - cross-references: rs1079598, rs1079597
    Observed alleles
        - A,C
        - A,T
        - G,C
        - G,T


--[ Marker Sequence ]--
>mh11KK-090
AGATTCGCCTTTC


--[ Target Amplicon Sequence with Alleles ]--
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
    ('mh04KK-010', 'mh04KK-010  MHDBM-07c8d144    GRCh38  chr4  1985210,1985244  2.7234  ALFRED'),
    ('mh08PK-46625', 'mh08PK-46625  MHDBM-840756f3    GRCh38  chr8  1194352,1194356,1194364,1194371  2.6273  10.1016/j.fsigen.2018.05.008'),
    ('mh04AT-10', 'mh04AT-10  MHDBM-07c8d144    GRCh38  chr4  1985210,1985244  3.8364  ISFG2019:P597'),
    ('mh01NH-03', 'mh01NH-03  MHDBM-e7a95c5e    GRCh38  chr1  184807944,184807966,184808042  2.3711  10.1016/j.legalmed.2015.06.003'),
    ('mh04CP-004', 'mh04CP-004  MHDBM-8408d717    GRCh38  chr4  7402842,7402854,7402870  2.992  10.1016/j.fsigen.2019.02.018')
])
def test_all_sources(name, data, capsys):
    marker = microhapdb.markers[microhapdb.markers.Name == name]
    microhapdb.marker.print_table(marker)
    terminal = capsys.readouterr()
    print(terminal.out)
    assert data in terminal.out
