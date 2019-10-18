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
    assert len(microhapdb.frequencies) == 81918 + 366 + 386 + 33


def test_allele_frequencies():
    """
    >>> import microhapdb
    >>> f = microhapdb.frequencies
    >>> f[f.Marker == 'mh15CP-003'].Allele.unique()
    array(['A,G,A', 'A,G,C', 'A,A,C', 'C,G,C'], dtype=object)
    >>> f[(f.Marker == 'mh15CP-003') & (f.Allele == 'A,A,C')]
               Marker Population Allele  Frequency
    60293  mh15CP-003  SA004046O  A,A,C      0.056
    60297  mh15CP-003  SA004047P  A,A,C      0.033
    60301  mh15CP-003  SA004048Q  A,A,C      0.028
    60305  mh15CP-003  SA004049R  A,A,C      0.318
    60309  mh15CP-003  SA004050J  A,A,C      0.286
    60313  mh15CP-003  SA004057Q  A,A,C      0.196
    60317  mh15CP-003  SA004058R  A,A,C      0.243
    60321  mh15CP-003  SA004059S  A,A,C      0.262
    60325  mh15CP-003  SA004060K  A,A,C      0.226
    60329  mh15CP-003  SA004108N  A,A,C      0.215
    60333  mh15CP-003  SA004109O  A,A,C      0.213
    60337  mh15CP-003  SA004110G  A,A,C      0.328
    60341  mh15CP-003  SA004111H  A,A,C      0.221
    60345  mh15CP-003  SA004238R  A,A,C      0.280
    60349  mh15CP-003  SA004239S  A,A,C      0.238
    60353  mh15CP-003  SA004240K  A,A,C      0.344
    60357  mh15CP-003  SA004241L  A,A,C      0.157
    60361  mh15CP-003  SA004242M  A,A,C      0.057
    60365  mh15CP-003  SA004243N  A,A,C      0.018
    60369  mh15CP-003  SA004244O  A,A,C      0.047
    60373  mh15CP-003  SA004245P  A,A,C      0.235
    60377  mh15CP-003  SA004246Q  A,A,C      0.272
    60381  mh15CP-003  SA004247R  A,A,C      0.245
    60385  mh15CP-003  SA004248S  A,A,C      0.040
    60389  mh15CP-003  SA004249T  A,A,C      0.268
    60393  mh15CP-003  SA004250L  A,A,C      0.293
    >>> f.query('Marker == "mh15CP-003" and Allele == "A,A,C" and Population == "SA004049R"')
               Marker Population Allele  Frequency
    60305  mh15CP-003  SA004049R  A,A,C      0.318
    """
    af = microhapdb.frequencies
    assert af.shape == (82703, 4)
    result = af[af.Marker == 'mh21KK-315'].Allele.unique()
    assert len(result) == 8
    result = af[(af.Marker == 'mh21KK-315') & (af.Allele == 'A,C,T')]
    assert len(result) == 96
    result = af.query('Marker == "mh21KK-315" & Allele == "A,C,T" & Population == "SA001773S"').Frequency.values[0]
    assert result == pytest.approx(0.025)
