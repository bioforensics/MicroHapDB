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
    assert len(microhapdb.frequencies) == 54347 + 366 + 137 + 33 + 103 + 58490


def test_allele_frequencies():
    """
    >>> import microhapdb
    >>> f = microhapdb.frequencies
    >>> f[f.Marker == 'mh15CP-003'].Allele.unique()
    array(['A,A,C', 'A,G,A', 'A,G,C', 'C,G,C'], dtype=object)
    >>> f[(f.Marker == 'mh15CP-003') & (f.Allele == 'A,A,C')]
               Marker Population Allele  Frequency
    84602  mh15CP-003        ACB  A,A,C      0.057
    84606  mh15CP-003        ASW  A,A,C      0.033
    84610  mh15CP-003        BEB  A,A,C      0.238
    84613  mh15CP-003        CDX  A,A,C      0.280
    84617  mh15CP-003        CEU  A,A,C      0.293
    84620  mh15CP-003        CHB  A,A,C      0.243
    84624  mh15CP-003        CHS  A,A,C      0.262
    84628  mh15CP-003        CLM  A,A,C      0.213
    84632  mh15CP-003        ESN  A,A,C      0.040
    84635  mh15CP-003        FIN  A,A,C      0.318
    84639  mh15CP-003        GBR  A,A,C      0.286
    84642  mh15CP-003        GIH  A,A,C      0.272
    84645  mh15CP-003        GWD  A,A,C      0.018
    84648  mh15CP-003        IBS  A,A,C      0.215
    84651  mh15CP-003        ITU  A,A,C      0.245
    84654  mh15CP-003        JPT  A,A,C      0.226
    84658  mh15CP-003        KHV  A,A,C      0.268
    84662  mh15CP-003        LWK  A,A,C      0.056
    84665  mh15CP-003        MSL  A,A,C      0.047
    84668  mh15CP-003        MXL  A,A,C      0.328
    84672  mh15CP-003        PEL  A,A,C      0.235
    84676  mh15CP-003        PJL  A,A,C      0.344
    84679  mh15CP-003        PUR  A,A,C      0.221
    84683  mh15CP-003        STU  A,A,C      0.157
    84687  mh15CP-003        TSI  A,A,C      0.196
    84690  mh15CP-003        YRI  A,A,C      0.028
    >>> f.query('Marker == "mh15CP-003" and Allele == "A,A,C" and Population == "FIN"')
               Marker Population Allele  Frequency
    84635  mh15CP-003        FIN  A,A,C      0.318
    """
    af = microhapdb.frequencies
    assert af.shape == (113476, 4)
    result = af[af.Marker == 'mh21KK-315'].Allele.unique()
    assert len(result) == 8
    result = af[(af.Marker == 'mh21KK-315') & (af.Allele == 'A,C,T')]
    assert len(result) == 81
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
