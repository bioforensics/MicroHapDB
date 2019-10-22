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
    assert len(microhapdb.frequencies) == 81918 + 366 + 386 + 33 + 103


def test_allele_frequencies():
    """
    >>> import microhapdb
    >>> f = microhapdb.frequencies
    >>> f[f.Marker == 'mh15CP-003'].Allele.unique()
    array(['A,G,A', 'A,G,C', 'A,A,C', 'C,G,C'], dtype=object)
    >>> f[(f.Marker == 'mh15CP-003') & (f.Allele == 'A,A,C')]
               Marker Population Allele  Frequency
    60373  mh15CP-003  SA004046O  A,A,C      0.056
    60377  mh15CP-003  SA004047P  A,A,C      0.033
    60381  mh15CP-003  SA004048Q  A,A,C      0.028
    60385  mh15CP-003  SA004049R  A,A,C      0.318
    60389  mh15CP-003  SA004050J  A,A,C      0.286
    60393  mh15CP-003  SA004057Q  A,A,C      0.196
    60397  mh15CP-003  SA004058R  A,A,C      0.243
    60401  mh15CP-003  SA004059S  A,A,C      0.262
    60405  mh15CP-003  SA004060K  A,A,C      0.226
    60409  mh15CP-003  SA004108N  A,A,C      0.215
    60413  mh15CP-003  SA004109O  A,A,C      0.213
    60417  mh15CP-003  SA004110G  A,A,C      0.328
    60421  mh15CP-003  SA004111H  A,A,C      0.221
    60425  mh15CP-003  SA004238R  A,A,C      0.280
    60429  mh15CP-003  SA004239S  A,A,C      0.238
    60433  mh15CP-003  SA004240K  A,A,C      0.344
    60437  mh15CP-003  SA004241L  A,A,C      0.157
    60441  mh15CP-003  SA004242M  A,A,C      0.057
    60445  mh15CP-003  SA004243N  A,A,C      0.018
    60449  mh15CP-003  SA004244O  A,A,C      0.047
    60453  mh15CP-003  SA004245P  A,A,C      0.235
    60457  mh15CP-003  SA004246Q  A,A,C      0.272
    60461  mh15CP-003  SA004247R  A,A,C      0.245
    60465  mh15CP-003  SA004248S  A,A,C      0.040
    60469  mh15CP-003  SA004249T  A,A,C      0.268
    60473  mh15CP-003  SA004250L  A,A,C      0.293
    >>> f.query('Marker == "mh15CP-003" and Allele == "A,A,C" and Population == "SA004049R"')
               Marker Population Allele  Frequency
    60385  mh15CP-003  SA004049R  A,A,C      0.318
    """
    af = microhapdb.frequencies
    assert af.shape == (82806, 4)
    result = af[af.Marker == 'mh21KK-315'].Allele.unique()
    assert len(result) == 8
    result = af[(af.Marker == 'mh21KK-315') & (af.Allele == 'A,C,T')]
    assert len(result) == 96
    result = af.query('Marker == "mh21KK-315" & Allele == "A,C,T" & Population == "SA001773S"').Frequency.values[0]
    assert result == pytest.approx(0.025)


@pytest.mark.parametrize('marker,pop,allele,data', [
    ('mh22KK-064', 'SA000009J', 'A,A,T,AATAATT', 'mh22KK-064  SA000009J  A,A,T,AATAATT      0.828'),
    ('mh06PK-24844', 'MHDBP-383d86606a', 'C,C,G,C,C,C,A,A,A,A', 'mh06PK-24844  MHDBP-383d86606a  C,C,G,C,C,C,A,A,A,A      0.005'),
    ('mh20AT-40', 'MHDBP-7c055e7ee8', 'T,C,G', 'mh20AT-40  MHDBP-7c055e7ee8  T,C,G     0.0806'),
    ('mh11NH-17', 'MHDBP-63967b883e', 'C,G,G', 'mh11NH-17  MHDBP-63967b883e  C,G,G      0.153'),
])
def test_all_sources(marker, pop, allele, data, capsys):
    arglist = ['frequency', '--marker', marker, '--population', pop, '--allele', allele]
    args = microhapdb.cli.get_parser().parse_args(arglist)
    microhapdb.cli.frequency.main(args)
    terminal = capsys.readouterr()
    print(terminal.out)
    assert data in terminal.out
