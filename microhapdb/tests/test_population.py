# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


import microhapdb
from microhapdb.population import standardize_ids
import pytest


def test_standardize_ids():
    assert list(standardize_ids(['SA004057Q']).values) == ['SA004057Q']
    assert list(standardize_ids(['Mende']).values) == ['SA004244O']
    assert list(standardize_ids(['Maya, Yucatan', 'SA000055K', 'Greeks']).values) == ['SA002767W', 'SA000055K', 'SA000013E']
    print(list(standardize_ids(['Han']).values))
    assert list(standardize_ids(['Han']).values) == ['SA004058R', 'SA004059S', 'MHDBP-48c2cfb2aa', 'SA000009J', 'SA000001B']


def test_assumptions():
    assert len(microhapdb.populations) == 96 + 3 + 1 + 1 + 1


def test_populations():
    """
    >>> import microhapdb
    >>> p = microhapdb.populations
    >>> p[p.ID == 'SA000040E']
               ID     Name  Source
    49  SA000040E  Kachari  ALFRED
    >>> p[p.ID == 'SA000936S']
               ID     Name  Source
    54  SA000936S  Koreans  ALFRED
    >>> p[p.Name == 'Han']
                      ID Name                        Source
    29         SA004058R  Han                        ALFRED
    30         SA004059S  Han                        ALFRED
    31  MHDBP-48c2cfb2aa  Han  10.1016/j.fsigen.2019.02.018
    32         SA000009J  Han                        ALFRED
    33         SA000001B  Han                        ALFRED
    >>> p.query('Name.str.contains("Afr")')
                     ID               Name                        Source
    1  MHDBP-3dab7bdd14             Africa  10.1016/j.fsigen.2018.05.008
    2         SA000101C  African Americans                        ALFRED
    3         SA004047P  African Americans                        ALFRED
    4         SA004242M    Afro-Caribbeans                        ALFRED
    """
    pop = microhapdb.populations
    assert pop.shape == (102, 3)
    assert pop[pop.ID == 'SA004049R'].Name.values == ['Finns']
    assert pop[pop.ID == 'SA000028K'].Name.values == ['Karitiana']
    result = pop[pop.Name.str.contains('Jews')].ID.values
    assert list(result) == ['SA000490N', 'SA000015G', 'SA000096P', 'SA000016H']


def test_pop_table(capsys):
    masai = microhapdb.populations[microhapdb.populations.Name == 'Masai']
    microhapdb.population.print_table(masai)
    testout = '''
        ID   Name  Source
 SA000854R  Masai  ALFRED
'''
    terminal = capsys.readouterr()
    assert terminal.out.strip() == testout.strip()


def test_pop_table_multi(capsys):
    hanchinese = microhapdb.populations[microhapdb.populations.Name == 'Han']
    microhapdb.population.print_table(hanchinese)
    testout = '''
               ID Name                        Source
        SA004058R  Han                        ALFRED
        SA004059S  Han                        ALFRED
 MHDBP-48c2cfb2aa  Han  10.1016/j.fsigen.2019.02.018
        SA000009J  Han                        ALFRED
        SA000001B  Han                        ALFRED
'''
    terminal = capsys.readouterr()
    assert terminal.out.strip() == testout.strip()


def test_pop_detail(capsys):
    hausa = microhapdb.populations[microhapdb.populations.Name == 'Hausa']
    microhapdb.population.print_detail(hausa)
    testout = '''
----------------------------------------------------------[ MicroHapulator ]----
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
'''
    terminal = capsys.readouterr()
    assert terminal.out.strip() == testout.strip()


def test_pop_detail_multi(capsys):
    japanese = microhapdb.populations[microhapdb.populations.Name == 'Japanese']
    microhapdb.population.print_detail(japanese)
    testout = '''
----------------------------------------------------------[ MicroHapulator ]----
Japanese    (SA004060K; source=ALFRED)

- 1070 total allele frequencies available
  for 198 markers

# Alleles | # Markers
---------------------
        19|*
        18|*
        14|**
        13|*
        12|******
        11|*
        10|*
         9|****
         8|**************
         7|********
         6|****************
         5|***********************************
         4|********************************************************************************************************
         2|****
--------------------------------------------------------------------------------

----------------------------------------------------------[ MicroHapulator ]----
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

----------------------------------------------------------[ MicroHapulator ]----
Japanese    (MHDBP-63967b883e; source=10.1016/j.legalmed.2015.06.003)

- 33 total allele frequencies available
  for 7 markers

# Alleles | # Markers
---------------------
         7|*
         6|*
         4|*****
--------------------------------------------------------------------------------
'''
    terminal = capsys.readouterr()
    assert terminal.out.strip() == testout.strip()


@pytest.mark.parametrize('ident,data', [
    ('SA000019K', 'SA000019K  Russians  ALFRED'),
    ('MHDBP-936bc36f79', 'MHDBP-936bc36f79  Asia  10.1016/j.fsigen.2018.05.008'),
    ('MHDBP-7c055e7ee8', 'MHDBP-7c055e7ee8  Swedish  ISFG2019:P597'),
    ('MHDBP-63967b883e', 'MHDBP-63967b883e  Japanese  10.1016/j.legalmed.2015.06.003'),
    ('MHDBP-48c2cfb2aa', 'MHDBP-48c2cfb2aa  Han  10.1016/j.fsigen.2019.02.018'),
])
def test_all_sources(ident, data, capsys):
    pop = microhapdb.populations[microhapdb.populations.ID == ident]
    microhapdb.population.print_table(pop)
    terminal = capsys.readouterr()
    print(terminal.out)
    assert data in terminal.out
