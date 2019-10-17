#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


import microhapdb
import pytest


def test_assumptions():
    assert len(microhapdb.populations) == 96 + 3 + 1
    assert len(microhapdb.markers) == 198 + 15 + 40


def test_allele_frequencies():
    """Allele frequencies for 198 markers across 96 populations.

    >>> import microhapdb
    >>> f = microhapdb.frequencies
    >>> f[f.Marker == 'mh15CP-003'].Allele.unique()
    array(['A,G,A', 'A,G,C', 'A,A,C', 'C,G,C'], dtype=object)
    >>> f[(f.Marker == 'mh15CP-003') & (f.Allele == 'A,A,C')]
               Marker Population Allele  Frequency
    60264  mh15CP-003  SA004046O  A,A,C      0.056
    60268  mh15CP-003  SA004047P  A,A,C      0.033
    60272  mh15CP-003  SA004048Q  A,A,C      0.028
    60276  mh15CP-003  SA004049R  A,A,C      0.318
    60280  mh15CP-003  SA004050J  A,A,C      0.286
    60284  mh15CP-003  SA004057Q  A,A,C      0.196
    60288  mh15CP-003  SA004058R  A,A,C      0.243
    60292  mh15CP-003  SA004059S  A,A,C      0.262
    60296  mh15CP-003  SA004060K  A,A,C      0.226
    60300  mh15CP-003  SA004108N  A,A,C      0.215
    60304  mh15CP-003  SA004109O  A,A,C      0.213
    60308  mh15CP-003  SA004110G  A,A,C      0.328
    60312  mh15CP-003  SA004111H  A,A,C      0.221
    60316  mh15CP-003  SA004238R  A,A,C      0.280
    60320  mh15CP-003  SA004239S  A,A,C      0.238
    60324  mh15CP-003  SA004240K  A,A,C      0.344
    60328  mh15CP-003  SA004241L  A,A,C      0.157
    60332  mh15CP-003  SA004242M  A,A,C      0.057
    60336  mh15CP-003  SA004243N  A,A,C      0.018
    60340  mh15CP-003  SA004244O  A,A,C      0.047
    60344  mh15CP-003  SA004245P  A,A,C      0.235
    60348  mh15CP-003  SA004246Q  A,A,C      0.272
    60352  mh15CP-003  SA004247R  A,A,C      0.245
    60356  mh15CP-003  SA004248S  A,A,C      0.040
    60360  mh15CP-003  SA004249T  A,A,C      0.268
    60364  mh15CP-003  SA004250L  A,A,C      0.293
    >>> f.query('Marker == "mh15CP-003" and Allele == "A,A,C" and Population == "SA004049R"')
               Marker Population Allele  Frequency
    60276  mh15CP-003  SA004049R  A,A,C      0.318
    """
    af = microhapdb.frequencies
    assert af.shape == (82670, 4)
    result = af[af.Marker == 'mh21KK-315'].Allele.unique()
    assert len(result) == 8
    result = af[(af.Marker == 'mh21KK-315') & (af.Allele == 'A,C,T')]
    assert len(result) == 96
    result = af.query('Marker == "mh21KK-315" & Allele == "A,C,T" & Population == "SA001773S"').Frequency.values[0]
    assert result == pytest.approx(0.025)


def test_markers():
    """Microhaplotype marker data

    >>> import microhapdb
    >>> m = microhapdb.markers
    >>> m[m.Name == 'mh18CP-005']
               Name          PermID Reference  Chrom                          Offsets   AvgAe  Source
    218  mh18CP-005  MHDBM-a85754d3    GRCh38  chr18  8892864,8892893,8892896,8892907  3.6614  ALFRED
    >>> m[m.Name == 'mh01KK-117']
              Name          PermID Reference Chrom                                  Offsets  AvgAe  Source
    14  mh01KK-117  MHDBM-39dc025f    GRCh38  chr1  204664211,204664268,204664371,204664397  3.933  ALFRED
    >>> m[m.Name == 'mh11PK-63643']
                 Name          PermID Reference  Chrom                                            Offsets  AvgAe                        Source
    135  mh11PK-63643  MHDBM-c5ce121f    GRCh38  chr11  34415814,34415816,34415818,34415835,34415836,3...  4.033  10.1016/j.fsigen.2018.05.008
    >>> m[m.Name == 'mh02AT-05']
             Name          PermID Reference Chrom                        Offsets   AvgAe         Source
    34  mh02AT-05  MHDBM-c3feaba8    GRCh38  chr2  160222899,160222923,160222938  5.1944  ISFG2019:P597
    """
    m = microhapdb.markers
    vm = microhapdb.variantmap
    assert m.shape == (253, 7)
    result = m[m.Chrom == 'chr19']
    assert len(result) == 6
    varids = vm[vm.Marker.isin(result.Name)].Variant.unique()
    assert len(varids) == 17


def test_populations():
    """Population data

    >>> import microhapdb
    >>> p = microhapdb.populations
    >>> p[p.ID == 'SA000040E']
               ID     Name  Source
    47  SA000040E  Kachari  ALFRED
    >>> p[p.ID == 'SA000936S']
               ID     Name  Source
    52  SA000936S  Koreans  ALFRED
    >>> p[p.Name == 'Han']
               ID Name  Source
    29  SA004059S  Han  ALFRED
    30  SA004058R  Han  ALFRED
    31  SA000001B  Han  ALFRED
    32  SA000009J  Han  ALFRED
    >>> p.query('Name.str.contains("Afr")')
              ID               Name                        Source
    1     Africa             Africa  10.1016/j.fsigen.2018.05.008
    2  SA004047P  African Americans                        ALFRED
    3  SA000101C  African Americans                        ALFRED
    4  SA004242M    Afro-Caribbeans                        ALFRED
    """
    pop = microhapdb.populations
    assert pop.shape == (100, 3)
    assert pop[pop.ID == 'SA004049R'].Name.values == ['Finns']
    assert pop[pop.ID == 'SA000028K'].Name.values == ['Karitiana']
    result = pop[pop.Name.str.contains('Jews')].ID.values
    assert list(result) == ['SA000490N', 'SA000015G', 'SA000096P', 'SA000016H']
