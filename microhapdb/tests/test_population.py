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
from microhapdb import Population
import pytest


def test_standardize_ids():
    assert Population.standardize_ids(["SA004057Q"]) == ["TSI"]
    assert Population.standardize_ids(["MSL"]) == ["MSL"]
    assert Population.standardize_ids(["Han"]) == ["MHDBP-48c2cfb2aa", "SA000001B", "SA000009J"]
    assert Population.standardize_ids(["Maya, Yucatan", "SA000055K", "Greeks"]) == [
        "SA000013E",
        "SA000055K",
        "SA002767W",
    ]


def test_assumptions():
    num_populations_per_source = [
        26,  # 1KGP
        70,  # ALFRED
        1,  # 10.1016/j.fsigen.2019.02.018
        7,  # 10.1016/j.fsigen.2020.102275
        1,  # 10.1016/j.legalmed.2015.06.003
        1,  # ISFG2019:P597
        3,  # 10.1016/j.fsigen.2018.05.008
    ]
    assert len(microhapdb.populations) == sum(num_populations_per_source)


def test_populations():
    """
    >>> from microhapdb import Population
    >>> pop = Population.from_id("SA000040E")
    >>> print(pop.popid, pop.name, pop.source)
    SA000040E Kachari ALFRED
    >>> Population.table_from_ids(["CEU", "IBS"])
          ID                                               Name Source
    40   IBS                        Iberian Population in Spain   1KGP
    103  CEU  Utah Residents (CEPH) with Northern and Wester...   1KGP
    >>> for pop in Population.from_query("Name.str.contains('Han')"):
    ...   print(pop.popid, pop.name, pop.source)
    MHDBP-48c2cfb2aa Han 10.1016/j.fsigen.2019.02.018
    SA000001B Han ALFRED
    SA000009J Han ALFRED
    CHB Han Chinese in Beijing, China 1KGP
    CHS Southern Han Chinese 1KGP
    >>> Population.table_from_query("Name.str.contains('Afr')")
                     ID                                     Name                        Source
    2  MHDBP-3dab7bdd14                                   Africa  10.1016/j.fsigen.2018.05.008
    3         SA000101C                        African Americans                        ALFRED
    4               ACB           African Caribbeans in Barbados                          1KGP
    5               ASW  Americans of African Ancestry in SW USA                          1KGP
    """
    pop = microhapdb.populations
    assert pop.shape == (109, 3)
    assert Population.from_id("FIN").name == "Finnish in Finland"
    assert Population.from_id("SA000028K").name == "Karitiana"
    result = Population.table_from_query("Name.str.contains('Jews')")
    assert result.ID.tolist() == ["SA000490N", "SA000015G", "SA000096P", "SA000016H"]


def test_pop_table():
    result = Population.table_from_ids(["Masai"])
    assert len(result) == 1
    assert result.ID.iloc[0] == "SA000854R"
    assert result.Source.iloc[0] == "ALFRED"


def test_pop_table_multi():
    result = Population.table_from_query("Name == 'Han'")
    assert len(result) == 3
    assert result.ID.tolist() == ["MHDBP-48c2cfb2aa", "SA000001B", "SA000009J"]


def test_pop_detail():
    pop = Population.from_id("Hausa")
    observed = pop.detail
    expected = """
--------------------------------------------------------------[ MicroHapDB ]----
Hausa    (SA000100B; source=ALFRED)

- 878 total allele frequencies available
  for 165 markers

# Alleles | # Markers
---------------------
        19|*
        18|*
        14|**
        13|*
        12|*****
        11|*
        10|*
         9|**
         8|************
         7|****
         6|************
         5|*********************
         4|**************************************************************************************************
         2|****
--------------------------------------------------------------------------------
"""
    assert observed.strip() == expected.strip()


def test_pop_detail_multi(capsys):
    for pop in Population.from_query("Name == 'Japanese'"):
        print(pop.detail)
    observed = capsys.readouterr().out
    expected = """
--------------------------------------------------------------[ MicroHapDB ]----
Japanese    (MHDBP-63967b883e; source=10.1016/j.legalmed.2015.06.003)

- 33 total allele frequencies available
  for 7 markers

# Alleles | # Markers
---------------------
         7|*
         6|*
         4|*****
--------------------------------------------------------------------------------

--------------------------------------------------------------[ MicroHapDB ]----
Japanese    (SA000010B; source=ALFRED)

- 878 total allele frequencies available
  for 165 markers

# Alleles | # Markers
---------------------
        19|*
        18|*
        14|**
        13|*
        12|*****
        11|*
        10|*
         9|**
         8|************
         7|****
         6|************
         5|*********************
         4|**************************************************************************************************
         2|****
--------------------------------------------------------------------------------
"""
    assert observed.strip() == expected.strip()


@pytest.mark.parametrize(
    "popid,name,source",
    [
        ("GBR", "British in England and Scotland", "1KGP"),
        ("SA000019K", "Russians", "ALFRED"),
        ("MHDBP-48c2cfb2aa", "Han", "10.1016/j.fsigen.2019.02.018"),
        ("mMHseq-Zaramo", "Zaramo", "10.1016/j.fsigen.2020.102275"),
        ("MHDBP-63967b883e", "Japanese", "10.1016/j.legalmed.2015.06.003"),
        ("MHDBP-7c055e7ee8", "Swedish", "ISFG2019:P597"),
        ("MHDBP-936bc36f79", "Asia", "10.1016/j.fsigen.2018.05.008"),
    ],
)
def test_all_sources(popid, name, source):
    pop = Population.from_id(popid)
    assert pop.popid == popid
    assert pop.name == name
    assert pop.source == source


def test_from_id_pop_not_found():
    with pytest.raises(ValueError, match=r"population 'Romulans' not found"):
        Population.from_id("Romulans")
