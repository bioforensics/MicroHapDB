# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


import microhapdb
from microhapdb.population import standardize_ids
import pytest


def test_assumptions():
    assert len(microhapdb.frequencies) == 81918 + 366 + 386


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
