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
        35041,  # byrska-bishop2022
        82167,  # kidd2018
        103,  # chen2019
        4737,  # gandotra2020
        33,  # hiroaki2015
        186,  # staadig2021
        366,  # vandergaag2018
    ]
    assert len(microhapdb.frequencies) == sum(num_allele_freqs_per_source)


def test_allele_frequencies():
    """
    >>> import microhapdb
    >>> f = microhapdb.frequencies
    >>> f[f.Marker == 'mh15CP-003'].Allele.unique()
    array(['A|A|C', 'A|G|A', 'A|G|C', 'C|G|C'], dtype=object)
    >>> f[(f.Marker == "mh15CP-003") & (f.Allele == "A|A|C")]
               Marker Population Allele  Frequency
    90512  mh15CP-003       1KGP  A|A|C    0.19947
    90518  mh15CP-003        ACB  A|A|C    0.05700
    90520  mh15CP-003        AFR  A|A|C    0.03565
    90525  mh15CP-003        ASW  A|A|C    0.03300
    90529  mh15CP-003        BEB  A|A|C    0.23800
    90533  mh15CP-003        CDX  A|A|C    0.28000
    90537  mh15CP-003        CEU  A|A|C    0.29300
    90541  mh15CP-003        CHB  A|A|C    0.24300
    90545  mh15CP-003        CHS  A|A|C    0.26200
    90549  mh15CP-003        CLM  A|A|C    0.21300
    90551  mh15CP-003        EAS  A|A|C    0.25828
    90557  mh15CP-003        ESN  A|A|C    0.04000
    90559  mh15CP-003        EUR  A|A|C    0.26046
    90565  mh15CP-003        FIN  A|A|C    0.31800
    90569  mh15CP-003        GBR  A|A|C    0.28600
    90573  mh15CP-003        GIH  A|A|C    0.27200
    90577  mh15CP-003        GWD  A|A|C    0.01800
    90581  mh15CP-003        IBS  A|A|C    0.21500
    90585  mh15CP-003        ITU  A|A|C    0.24500
    90589  mh15CP-003        JPT  A|A|C    0.22600
    90593  mh15CP-003        KHV  A|A|C    0.26800
    90597  mh15CP-003        LWK  A|A|C    0.05600
    90601  mh15CP-003        MSL  A|A|C    0.04700
    90605  mh15CP-003        MXL  A|A|C    0.32800
    90609  mh15CP-003        PEL  A|A|C    0.23500
    90613  mh15CP-003        PJL  A|A|C    0.34400
    90617  mh15CP-003        PUR  A|A|C    0.22100
    90619  mh15CP-003        SAS  A|A|C    0.24806
    90625  mh15CP-003        STU  A|A|C    0.15700
    90629  mh15CP-003        TSI  A|A|C    0.19600
    90633  mh15CP-003        YRI  A|A|C    0.02800
    >>> f.query("Marker == 'mh15CP-003' and Allele == 'A|A|C' and Population == 'FIN'")
               Marker Population Allele  Frequency
    90565  mh15CP-003        FIN  A|A|C      0.318
    """
    af = microhapdb.frequencies
    assert af.shape == (122633, 4)
    result = af[af.Marker == "mh21KK-315.v1"].Allele.unique()
    assert len(result) == 8
    result = af[(af.Marker == "mh21KK-315.v1") & (af.Allele == "A|C|T")]
    assert len(result) == 101
    result = af.query(
        'Marker == "mh21KK-315.v1" & Allele == "A|C|T" & Population == "SA001773S"'
    ).Frequency.values[0]
    assert result == pytest.approx(0.025)


@pytest.mark.parametrize(
    "marker,pop,allele,frequency",
    [
        ("mh10USC-10qC", "EUR", "G|C|A", 0.02567),
        (
            "mh11KK-089",
            "SA000009J",
            "A|T",
            0.482,
        ),
        ("mh01CP-016", "MHDBP-48c2cfb2aa", "T|G|A", 0.2916),
        ("mh18KK-293.v2", "mMHseq-Chagga", "G|C|G|A|T|A|G", 0.011),
        ("mh11NH-17", "MHDBP-63967b883e", "C|G|G", 0.153),
        ("mh20KK-307.v3", "MHDBP-7c055e7ee8", "T|C|G", 0.127),
        (
            "mh06PK-24844",
            "MHDBP-383d86606a",
            "C|C|G|C|C|C|A|A|A|A",
            0.005,
        ),
    ],
)
def test_all_sources(marker, pop, allele, frequency):
    freq = microhapdb.frequencies
    result = freq[(freq.Marker == marker) & (freq.Population == pop) & (freq.Allele == allele)]
    assert len(result) == 1
    assert result.Frequency.iloc[0] == pytest.approx(frequency)
