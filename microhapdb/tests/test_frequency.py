# -------------------------------------------------------------------------------------------------
# Copyright (c) 2019, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------


import microhapdb
import pytest


def test_assumptions():
    num_allele_freqs_per_source = [
        89524,  # 1KGP
        54347,  # ALFRED
        103,  # 10.1016/j.fsigen.2019.02.018
        4737,  # 10.1016/j.fsigen.2020.102275
        33,  # 10.1016/j.legalmed.2015.06.003
        137,  # ISFG2019:P597
        366,  # 10.1016/j.fsigen.2018.05.008
    ]
    assert len(microhapdb.frequencies) == sum(num_allele_freqs_per_source)


def test_allele_frequencies():
    """
    >>> import microhapdb
    >>> f = microhapdb.frequencies
    >>> f[f.Marker == 'mh15CP-003'].Allele.unique()
    array(['A,A,C', 'A,G,A', 'A,G,C', 'C,G,C'], dtype=object)
    >>> f[(f.Marker == "mh15CP-003") & (f.Allele == "A,A,C")]
                Marker Population Allele  Frequency
    111073  mh15CP-003        ACB  A,A,C      0.057
    111077  mh15CP-003        ASW  A,A,C      0.033
    111081  mh15CP-003        BEB  A,A,C      0.238
    111084  mh15CP-003        CDX  A,A,C      0.280
    111088  mh15CP-003        CEU  A,A,C      0.293
    111091  mh15CP-003        CHB  A,A,C      0.243
    111095  mh15CP-003        CHS  A,A,C      0.262
    111099  mh15CP-003        CLM  A,A,C      0.213
    111103  mh15CP-003        ESN  A,A,C      0.040
    111106  mh15CP-003        FIN  A,A,C      0.318
    111110  mh15CP-003        GBR  A,A,C      0.286
    111113  mh15CP-003        GIH  A,A,C      0.272
    111116  mh15CP-003        GWD  A,A,C      0.018
    111119  mh15CP-003        IBS  A,A,C      0.215
    111122  mh15CP-003        ITU  A,A,C      0.245
    111125  mh15CP-003        JPT  A,A,C      0.226
    111129  mh15CP-003        KHV  A,A,C      0.268
    111133  mh15CP-003        LWK  A,A,C      0.056
    111136  mh15CP-003        MSL  A,A,C      0.047
    111139  mh15CP-003        MXL  A,A,C      0.328
    111143  mh15CP-003        PEL  A,A,C      0.235
    111147  mh15CP-003        PJL  A,A,C      0.344
    111150  mh15CP-003        PUR  A,A,C      0.221
    111154  mh15CP-003        STU  A,A,C      0.157
    111158  mh15CP-003        TSI  A,A,C      0.196
    111161  mh15CP-003        YRI  A,A,C      0.028
    >>> f.query("Marker == 'mh15CP-003' and Allele == 'A,A,C' and Population == 'FIN'")
                Marker Population Allele  Frequency
    111106  mh15CP-003        FIN  A,A,C      0.318
    """
    af = microhapdb.frequencies
    assert af.shape == (149247, 4)
    result = af[af.Marker == "mh21KK-315"].Allele.unique()
    assert len(result) == 8
    result = af[(af.Marker == "mh21KK-315") & (af.Allele == "A,C,T")]
    assert len(result) == 81
    result = af.query(
        'Marker == "mh21KK-315" & Allele == "A,C,T" & Population == "SA001773S"'
    ).Frequency.values[0]
    assert result == pytest.approx(0.025)


@pytest.mark.parametrize(
    "marker,pop,allele,frequency",
    [
        ("mh10USC-10qC", "MXL", "G,C,A", 0.039),  # 1KGP
        (
            "mh22KK-064",
            "SA000009J",
            "A,A,T,AATAATT",
            0.828,
        ),  # ALFRED
        ("mh01CP-016", "MHDBP-48c2cfb2aa", "T,G,A", 0.2916),  # 10.1016/j.fsigen.2019.02.018
        ("mh18KKCS-293", "mMHseq-Chagga", "G,C,G,A,T,A,G", 0.011),  # 10.1016/j.fsigen.2020.102275
        ("mh11NH-17", "MHDBP-63967b883e", "C,G,G", 0.153),  # 10.1016/j.legalmed.2015.06.003
        ("mh20AT-40", "MHDBP-7c055e7ee8", "T,C,G", 0.0806),  # ISFG2019:P597
        (
            "mh06PK-24844",
            "MHDBP-383d86606a",
            "C,C,G,C,C,C,A,A,A,A",
            0.005,
        ),  # 10.1016/j.fsigen.2018.05.008
    ],
)
def test_all_sources(marker, pop, allele, frequency):
    freq = microhapdb.frequencies
    result = freq[(freq.Marker == marker) & (freq.Population == pop) & (freq.Allele == allele)]
    assert len(result) == 1
    assert result.Frequency.iloc[0] == pytest.approx(frequency)
