# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


from io import StringIO
import microhapdb
import pandas
import pytest


def test_assumptions():
    assert len(microhapdb.frequencies) == 54347 + 366 + 137 + 33 + 103 + 66565 + 4737


def test_allele_frequencies():
    """
    >>> import microhapdb
    >>> f = microhapdb.frequencies
    >>> f[f.Marker == 'mh15CP-003'].Allele.unique()
    array(['A,A,C', 'A,G,A', 'A,G,C', 'C,G,C'], dtype=object)
    >>> f[(f.Marker == 'mh15CP-003') & (f.Allele == 'A,A,C')]
               Marker Population Allele  Frequency
    93520  mh15CP-003        ACB  A,A,C      0.057
    93524  mh15CP-003        ASW  A,A,C      0.033
    93528  mh15CP-003        BEB  A,A,C      0.238
    93531  mh15CP-003        CDX  A,A,C      0.280
    93535  mh15CP-003        CEU  A,A,C      0.293
    93538  mh15CP-003        CHB  A,A,C      0.243
    93542  mh15CP-003        CHS  A,A,C      0.262
    93546  mh15CP-003        CLM  A,A,C      0.213
    93550  mh15CP-003        ESN  A,A,C      0.040
    93553  mh15CP-003        FIN  A,A,C      0.318
    93557  mh15CP-003        GBR  A,A,C      0.286
    93560  mh15CP-003        GIH  A,A,C      0.272
    93563  mh15CP-003        GWD  A,A,C      0.018
    93566  mh15CP-003        IBS  A,A,C      0.215
    93569  mh15CP-003        ITU  A,A,C      0.245
    93572  mh15CP-003        JPT  A,A,C      0.226
    93576  mh15CP-003        KHV  A,A,C      0.268
    93580  mh15CP-003        LWK  A,A,C      0.056
    93583  mh15CP-003        MSL  A,A,C      0.047
    93586  mh15CP-003        MXL  A,A,C      0.328
    93590  mh15CP-003        PEL  A,A,C      0.235
    93594  mh15CP-003        PJL  A,A,C      0.344
    93597  mh15CP-003        PUR  A,A,C      0.221
    93601  mh15CP-003        STU  A,A,C      0.157
    93605  mh15CP-003        TSI  A,A,C      0.196
    93608  mh15CP-003        YRI  A,A,C      0.028
    >>> f.query('Marker == "mh15CP-003" and Allele == "A,A,C" and Population == "FIN"')
               Marker Population Allele  Frequency
    93553  mh15CP-003        FIN  A,A,C      0.318
    """
    af = microhapdb.frequencies
    assert af.shape == (126288, 4)
    result = af[af.Marker == 'mh21KK-315'].Allele.unique()
    assert len(result) == 8
    result = af[(af.Marker == 'mh21KK-315') & (af.Allele == 'A,C,T')]
    assert len(result) == 81
    result = af.query('Marker == "mh21KK-315" & Allele == "A,C,T" & Population == "SA001773S"').Frequency.values[0]
    assert result == pytest.approx(0.025)


@pytest.mark.parametrize('marker,pop,allele,data', [
    ('mh22KK-064', 'SA000009J', 'A,A,T,AATAATT', 'mh22KK-064  SA000009J A,A,T,AATAATT      0.828'),
    ('mh06PK-24844', 'MHDBP-383d86606a', 'C,C,G,C,C,C,A,A,A,A', 'mh06PK-24844 MHDBP-383d86606a C,C,G,C,C,C,A,A,A,A      0.005'),
    ('mh20AT-40', 'MHDBP-7c055e7ee8', 'T,C,G', 'mh20AT-40 MHDBP-7c055e7ee8  T,C,G     0.0806'),
    ('mh11NH-17', 'MHDBP-63967b883e', 'C,G,G', 'mh11NH-17 MHDBP-63967b883e  C,G,G      0.153'),
    ('mh01CP-016', 'MHDBP-48c2cfb2aa', 'T,G,A', 'mh01CP-016 MHDBP-48c2cfb2aa  T,G,A     0.2916'),
])
def test_all_sources(marker, pop, allele, data, capsys):
    arglist = ['frequency', '--marker', marker, '--population', pop, '--allele', allele]
    args = microhapdb.cli.get_parser().parse_args(arglist)
    microhapdb.cli.frequency.main(args)
    terminal = capsys.readouterr()
    print(terminal.out)
    assert data in terminal.out


def test_mhpl8r(capsys):
    arglist = ['frequency', '--marker', 'mh02USC-2pA', '--population', 'JPT', '--format', 'mhpl8r']
    args = microhapdb.cli.get_parser().parse_args(arglist)
    microhapdb.cli.frequency.main(args)
    terminal = capsys.readouterr()
    result = pandas.read_csv(StringIO(terminal.out), sep='\t')
    assert result.shape == (4, 3)
    assert result.Haplotype.iloc[0] == 'A,A,G,A'
    assert result.Frequency.iloc[0] == pytest.approx(0.005)


def test_mhpl8r_multi_pop(capsys):
    arglist = ['frequency', '--marker', 'mh02USC-2pA', '--format', 'mhpl8r']
    args = microhapdb.cli.get_parser().parse_args(arglist)
    microhapdb.cli.frequency.main(args)
    terminal = capsys.readouterr()
    assert 'warning: frequencies for 26 populations recovered, expected only 1' in terminal.err


def test_bad_format():
    arglist = ['frequency', '--marker', 'mh02USC-2pA', '--population', 'JPT', '--format', 'detail']
    args = microhapdb.cli.get_parser().parse_args(arglist)
    with pytest.raises(NotImplementedError):
        microhapdb.cli.frequency.main(args)
    args.format = 'BoGuS'
    with pytest.raises(ValueError, match=r'unsupported view format "BoGuS"'):
        microhapdb.cli.frequency.main(args)
