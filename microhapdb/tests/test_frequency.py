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
from random import choice


def test_assumptions():
    num_allele_freqs_per_source = [
        856067,  # Byrska-Bishop2022
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
                Marker Population Allele  Frequency  Count             Source
    716565  mh15CP-003        YRI  A|A|C    0.02893    242  Byrska-Bishop2022
    716566  mh15CP-003        YRI  A|G|A    0.59091    242  Byrska-Bishop2022
    716567  mh15CP-003        YRI  A|G|C    0.38017    242  Byrska-Bishop2022
    716568  mh15CP-003        YRI  A|G|A    0.57400    216           Kidd2018
    716569  mh15CP-003        YRI  A|G|C    0.39800    216           Kidd2018
    716570  mh15CP-003        YRI  A|A|C    0.02800    216           Kidd2018
    716571  mh15CP-003        YRI  C|G|C    0.00000    216           Kidd2018
    >>> f.query("Marker == 'mh15CP-003' and Allele == 'A|A|C' and Population == 'FIN'")
                Marker Population Allele  Frequency  Count             Source
    716442  mh15CP-003        FIN  A|A|C    0.31818    198  Byrska-Bishop2022
    716448  mh15CP-003        FIN  A|A|C    0.31800    198           Kidd2018
    """
    af = microhapdb.frequencies
    assert af.shape == (944736, 6)
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


def test_all_consistent():
    inconsistent_markers = set()
    for markerid, table in microhapdb.frequencies.groupby("Marker"):
        numvars = set([a.count("|") + 1 for a in table.Allele])
        if len(numvars) > 1:
            inconsistent_markers.add(markerid)
    assert len(inconsistent_markers) == 0, sorted(inconsistent_markers)


@pytest.mark.parametrize(
    "marker,population,counts",
    [
        ("mh01CP-012", "EUR", [1052] * 4),
        ("mh11HYP-29", "ASW", [122] * 3),
        ("mh03USC-3qB", "CHS", [224] * 4),
    ],
)
def test_counts_simple(marker, population, counts):
    freq = microhapdb.frequencies
    subset = freq[(freq.Marker == marker) & (freq.Population == population)]
    assert subset.Count.to_list() == counts


def test_counts_random():
    freq = microhapdb.frequencies
    pops = microhapdb.populations
    for _ in range(5):
        random_marker = choice(microhapdb.markers.Name)
        random_population = choice(pops[pops.Source == "Byrska-Bishop2022"].ID.to_list())
        subset = freq[
            (freq.Marker == random_marker)
            & (freq.Population == random_population)
            & (freq.Source == "Byrska-Bishop2022")
        ]
        assert len(set(subset.Count)) == 1, subset
