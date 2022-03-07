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
from tempfile import NamedTemporaryFile


def test_assumptions():
    assert len(microhapdb.frequencies) == 54347 + 366 + 137 + 33 + 103 + 66565 + 4737 + 4131


def test_allele_frequencies():
    """
    >>> import microhapdb
    >>> f = microhapdb.frequencies
    >>> f[f.Marker == 'mh15CP-003'].Allele.unique()
    array(['A,A,C', 'A,G,A', 'A,G,C', 'C,G,C'], dtype=object)
    >>> f[(f.Marker == 'mh15CP-003') & (f.Allele == 'A,A,C')]
               Marker Population Allele  Frequency
    96569  mh15CP-003        ACB  A,A,C      0.057
    96573  mh15CP-003        ASW  A,A,C      0.033
    96577  mh15CP-003        BEB  A,A,C      0.238
    96580  mh15CP-003        CDX  A,A,C      0.280
    96584  mh15CP-003        CEU  A,A,C      0.293
    96587  mh15CP-003        CHB  A,A,C      0.243
    96591  mh15CP-003        CHS  A,A,C      0.262
    96595  mh15CP-003        CLM  A,A,C      0.213
    96599  mh15CP-003        ESN  A,A,C      0.040
    96602  mh15CP-003        FIN  A,A,C      0.318
    96606  mh15CP-003        GBR  A,A,C      0.286
    96609  mh15CP-003        GIH  A,A,C      0.272
    96612  mh15CP-003        GWD  A,A,C      0.018
    96615  mh15CP-003        IBS  A,A,C      0.215
    96618  mh15CP-003        ITU  A,A,C      0.245
    96621  mh15CP-003        JPT  A,A,C      0.226
    96625  mh15CP-003        KHV  A,A,C      0.268
    96629  mh15CP-003        LWK  A,A,C      0.056
    96632  mh15CP-003        MSL  A,A,C      0.047
    96635  mh15CP-003        MXL  A,A,C      0.328
    96639  mh15CP-003        PEL  A,A,C      0.235
    96643  mh15CP-003        PJL  A,A,C      0.344
    96646  mh15CP-003        PUR  A,A,C      0.221
    96650  mh15CP-003        STU  A,A,C      0.157
    96654  mh15CP-003        TSI  A,A,C      0.196
    96657  mh15CP-003        YRI  A,A,C      0.028
    >>> f.query('Marker == "mh15CP-003" and Allele == "A,A,C" and Population == "FIN"')
               Marker Population Allele  Frequency
    96602  mh15CP-003        FIN  A,A,C      0.318
    """
    af = microhapdb.frequencies
    assert af.shape == (130419, 4)
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


def test_mhpl8r_panel(capsys):
    with NamedTemporaryFile(mode="wt") as tempfile:
        print('mh02USC-2pA', 'mh08USC-8qA', 'mh17USC-17qA', sep='\n', file=tempfile, flush=True)
        arglist = ['frequency', '--panel', tempfile.name, '--population', 'GBR', '--format', 'mhpl8r']
        args = microhapdb.cli.get_parser().parse_args(arglist)
        microhapdb.cli.frequency.main(args)
    terminal = capsys.readouterr()
    result = pandas.read_csv(StringIO(terminal.out), sep='\t')
    print(result)
    assert result.shape == (13, 3)
    assert result.Haplotype.iloc[7] == 'A,G,T'
    assert result.Frequency.iloc[7] == pytest.approx(0.429)


def test_mhpl8r_multi_pop(capsys):
    arglist = [
        'frequency', '--marker', 'mh02USC-2pA', 'mh08USC-8qA', 'mh17USC-17qA',
        '--population', 'CLM', 'GIH', 'ASW', 'CEU', '--format', 'mhpl8r'
    ]
    args = microhapdb.cli.get_parser().parse_args(arglist)
    microhapdb.cli.frequency.main(args)
    terminal = capsys.readouterr()
    assert 'warning: frequencies for 4 populations recovered, expected only 1' in terminal.err


def test_efm(capsys):
    arglist = [
        'frequency', '--population=CEU', '--format=efm', '--marker', 'mh01USC-1pD', 'mh17USC-17pA',
        'mh15USC-15qA'
    ]
    args = microhapdb.cli.get_parser().parse_args(arglist)
    microhapdb.cli.frequency.main(args)
    terminal = capsys.readouterr()
    result = pandas.read_csv(StringIO(terminal.out))
    print(result)
    assert result.shape == (13, 4)
    assert result['Allele'].iloc[3] == "C,T,C"
    assert result['mh01USC-1pD'].iloc[3] == pytest.approx(0.101)
    assert pandas.isna(result['mh15USC-15qA'].iloc[3])
    assert result['mh17USC-17pA'].iloc[3] == pytest.approx(0.02)


def test_efm_multi_pop():
    arglist = [
        'frequency', '--population', 'CEU', 'IBS', '--format=efm', '--marker', 'mh01USC-1pD',
        'mh17USC-17pA', 'mh15USC-15qA'
    ]
    args = microhapdb.cli.get_parser().parse_args(arglist)
    with pytest.raises(ValueError, match="must specify one and only one population with --format=efm"):
        microhapdb.cli.frequency.main(args)


def test_bad_format():
    arglist = ['frequency', '--marker', 'mh02USC-2pA', '--population', 'JPT', '--format', 'detail']
    args = microhapdb.cli.get_parser().parse_args(arglist)
    with pytest.raises(NotImplementedError):
        microhapdb.cli.frequency.main(args)
    args.format = 'BoGuS'
    with pytest.raises(ValueError, match=r'unsupported view format "BoGuS"'):
        microhapdb.cli.frequency.main(args)
