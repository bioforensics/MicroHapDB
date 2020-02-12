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
    array(['A,A,C', 'A,G,A', 'A,G,C', 'C,G,C'], dtype=object)
    >>> f[(f.Marker == 'mh15CP-003') & (f.Allele == 'A,A,C')]
               Marker Population Allele  Frequency
    84294  mh15CP-003        ACB  A,A,C      0.057
    84298  mh15CP-003        ASW  A,A,C      0.033
    84302  mh15CP-003        BEB  A,A,C      0.238
    84305  mh15CP-003        CDX  A,A,C      0.280
    84309  mh15CP-003        CEU  A,A,C      0.293
    84312  mh15CP-003        CHB  A,A,C      0.243
    84316  mh15CP-003        CHS  A,A,C      0.262
    84320  mh15CP-003        CLM  A,A,C      0.213
    84324  mh15CP-003        ESN  A,A,C      0.040
    84327  mh15CP-003        FIN  A,A,C      0.318
    84331  mh15CP-003        GBR  A,A,C      0.286
    84334  mh15CP-003        GIH  A,A,C      0.272
    84337  mh15CP-003        GWD  A,A,C      0.018
    84340  mh15CP-003        IBS  A,A,C      0.215
    84343  mh15CP-003        ITU  A,A,C      0.245
    84346  mh15CP-003        JPT  A,A,C      0.226
    84350  mh15CP-003        KHV  A,A,C      0.268
    84354  mh15CP-003        LWK  A,A,C      0.056
    84357  mh15CP-003        MSL  A,A,C      0.047
    84360  mh15CP-003        MXL  A,A,C      0.328
    84364  mh15CP-003        PEL  A,A,C      0.235
    84368  mh15CP-003        PJL  A,A,C      0.344
    84371  mh15CP-003        PUR  A,A,C      0.221
    84375  mh15CP-003        STU  A,A,C      0.157
    84379  mh15CP-003        TSI  A,A,C      0.196
    84382  mh15CP-003        YRI  A,A,C      0.028
    >>> f.query('Marker == "mh15CP-003" and Allele == "A,A,C" and Population == "FIN"')
               Marker Population Allele  Frequency
    84327  mh15CP-003        FIN  A,A,C      0.318
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
    ('mh01CP-016', 'MHDBP-48c2cfb2aa', 'T,G,A', 'mh01CP-016  MHDBP-48c2cfb2aa  T,G,A     0.2916'),
])
def test_all_sources(marker, pop, allele, data, capsys):
    arglist = ['frequency', '--marker', marker, '--population', pop, '--allele', allele]
    args = microhapdb.cli.get_parser().parse_args(arglist)
    microhapdb.cli.frequency.main(args)
    terminal = capsys.readouterr()
    print(terminal.out)
    assert data in terminal.out
