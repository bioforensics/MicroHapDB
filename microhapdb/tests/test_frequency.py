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
        89524,  # Auton2015
        35041,  # Byrska-Bishop2022
        103,  # Chen2019
        4737,  # Gandotra2020
        33,  # Hiroaki2015
        82167,  # Kidd2018
        186,  # Staadig2021
        366,  # vanderGaag2018
    ]
    assert len(microhapdb.frequencies) == sum(num_allele_freqs_per_source)


def test_allele_frequencies():
    """
    >>> import microhapdb
    >>> f = microhapdb.frequencies
    >>> f[f.Marker == 'mh15CP-003'].Allele.unique()
    array(['A|A|C', 'A|G|A', 'A|G|C', 'C|G|C'], dtype=object)
    >>> f[(f.Marker == "mh15CP-003") & (f.Allele == "A|A|C") & (f.Source == "Byrska-Bishop2022")]
                Marker Population Allele  Frequency             Source
    158523  mh15CP-003       1KGP  A|A|C    0.19947  Byrska-Bishop2022
    158535  mh15CP-003        AFR  A|A|C    0.03565  Byrska-Bishop2022
    158592  mh15CP-003        EAS  A|A|C    0.25828  Byrska-Bishop2022
    158603  mh15CP-003        EUR  A|A|C    0.26046  Byrska-Bishop2022
    158711  mh15CP-003        SAS  A|A|C    0.24806  Byrska-Bishop2022
    >>> f.query("Marker == 'mh15CP-003' and Allele == 'A|A|C' and Population == 'FIN'")
                Marker Population Allele  Frequency     Source
    158607  mh15CP-003        FIN  A|A|C      0.318  Auton2015
    158613  mh15CP-003        FIN  A|A|C      0.318   Kidd2018
    """
    af = microhapdb.frequencies
    assert af.shape == (212157, 5)
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
