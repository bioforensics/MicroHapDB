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
    assert Population.standardize_ids(["SA002766V"]) == ["SA002766V"]
    assert Population.standardize_ids(["Samoans"]) == ["SA000072J"]
    assert Population.standardize_ids(["Han"]) == ["MHDBP-48c2cfb2aa", "SA000001B", "SA000009J"]
    assert Population.standardize_ids(["Maya, Yucatan", "SA000055K", "Greeks"]) == [
        "SA000013E",
        "SA000055K",
        "SA002767W",
    ]


def test_assumptions():
    num_populations_per_source = [
        31,  # Byrska-Bishop2022
        1,  # Chen2019
        7,  # Gandotra2020
        1,  # Hiroaki2015
        70,  # Kidd2018
        1,  # Staadig2021
        1,  # Turchi2019
        10,  # Zou2022
        3,  # vanderGaag2018
    ]
    assert len(microhapdb.populations) == sum(num_populations_per_source)


def test_populations():
    """
    >>> pop = Population.from_id("SA000040E")
    >>> print(pop.popid, pop.name, pop.source)
    SA000040E Kachari Kidd2018
    >>> Population.table_from_ids(["EAS", "SAS"])
          ID        Name             Source
    27   EAS   East Asia  Byrska-Bishop2022
    104  SAS  South Asia  Byrska-Bishop2022
    >>> for pop in Population.from_query("Name.str.contains('Han')"):
    ...   print(pop.popid, pop.name, pop.source)
    ChengduHan Chengdu Han Zou2022
    HainanHan Hainan Han Zou2022
    MHDBP-48c2cfb2aa Han Chen2019
    SA000001B Han Kidd2018
    SA000009J Han Kidd2018
    CHB Han Chinese in Beijing, China Byrska-Bishop2022
    CHS Southern Han Chinese Byrska-Bishop2022
    >>> Population.table_from_query("Name.str.contains('Afr')")
                     ID                                     Name             Source
    3               AFR                                   Africa  Byrska-Bishop2022
    4  MHDBP-3dab7bdd14                                   Africa     vanderGaag2018
    5         SA000101C                        African Americans           Kidd2018
    6               ACB           African Caribbeans in Barbados  Byrska-Bishop2022
    7               ASW  Americans of African Ancestry in SW USA  Byrska-Bishop2022
    """
    pop = microhapdb.populations
    assert pop.shape == (125, 3)
    assert Population.from_id("MHDBP-7c055e7ee8").name == "Swedish"
    assert Population.from_id("SA000028K").name == "Karitiana"
    result = Population.table_from_query("Name.str.contains('Jews')")
    assert result.ID.tolist() == ["SA000490N", "SA000015G", "SA000096P", "SA000016H"]


def test_pop_table():
    result = Population.table_from_ids(["Masai"])
    assert len(result) == 1
    assert result.ID.iloc[0] == "SA000854R"
    assert result.Source.iloc[0] == "Kidd2018"


def test_pop_table_multi():
    result = Population.table_from_query("Name == 'Han'")
    assert len(result) == 3
    assert sorted(result.ID.tolist()) == ["MHDBP-48c2cfb2aa", "SA000001B", "SA000009J"]


def test_pop_detail():
    pop = Population.from_id("Hausa")
    observed = pop.detail
    expected = """
--------------------------------------------------------------[ MicroHapDB ]----
Hausa    (SA000100B; source=Kidd2018)

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
Japanese    (MHDBP-63967b883e; source=Hiroaki2015)

- 33 total allele frequencies available
  for 7 markers

# Alleles | # Markers
---------------------
         7|*
         6|*
         4|*****
--------------------------------------------------------------------------------

--------------------------------------------------------------[ MicroHapDB ]----
Japanese    (SA000010B; source=Kidd2018)

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
        ("SAS", "South Asia", "Byrska-Bishop2022"),
        ("SA000019K", "Russians", "Kidd2018"),
        ("MHDBP-48c2cfb2aa", "Han", "Chen2019"),
        ("mMHseq-Zaramo", "Zaramo", "Gandotra2020"),
        ("MHDBP-63967b883e", "Japanese", "Hiroaki2015"),
        ("MHDBP-7c055e7ee8", "Swedish", "Staadig2021"),
        ("MHDBP-936bc36f79", "Asia", "vanderGaag2018"),
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
