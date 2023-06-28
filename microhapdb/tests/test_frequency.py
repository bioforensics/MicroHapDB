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
        526413,  # Byrska-Bishop2022
        103,  # Chen2019
        4737,  # Gandotra2020
        33,  # Hiroaki2015
        82167,  # Kidd2018
        186,  # Staadig2021
        427,  # Turchi2019
        650,  # Zou2022
        366,  # vanderGaag2018
    ]
    assert len(microhapdb.frequencies) == sum(num_allele_freqs_per_source)


def test_allele_frequencies():
    """
    >>> f = microhapdb.frequencies
    >>> f[(f.Marker == "mh15CP-003") & (f.Population == "YRI")]
                Marker Population Allele  Frequency             Source
    474187  mh15CP-003        YRI  A|A|C    0.02893  Byrska-Bishop2022
    474188  mh15CP-003        YRI  A|G|A    0.59091  Byrska-Bishop2022
    474189  mh15CP-003        YRI  A|G|C    0.38017  Byrska-Bishop2022
    474190  mh15CP-003        YRI  A|G|A    0.57400           Kidd2018
    474191  mh15CP-003        YRI  A|G|C    0.39800           Kidd2018
    474192  mh15CP-003        YRI  A|A|C    0.02800           Kidd2018
    474193  mh15CP-003        YRI  C|G|C    0.00000           Kidd2018
    >>> f.query("Marker == 'mh15CP-003' and Allele == 'A|A|C' and Population == 'FIN'")
                Marker Population Allele  Frequency             Source
    474064  mh15CP-003        FIN  A|A|C    0.31818  Byrska-Bishop2022
    474070  mh15CP-003        FIN  A|A|C    0.31800           Kidd2018
    """
    af = microhapdb.frequencies
    assert af.shape == (615082, 5)
    result = af[af.Marker == "mh21KK-315.v1"].Allele.unique()
    assert len(result) == 8
    result = af[(af.Marker == "mh21KK-315.v1") & (af.Allele == "A|C|T")]
    assert len(result) == 115
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
