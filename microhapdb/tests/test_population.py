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
    assert list(standardize_ids(['Han']).values) == ['SA004059S', 'SA004058R', 'SA000001B', 'SA000009J']


def test_assumptions():
    assert len(microhapdb.populations) == 96 + 3 + 1 + 1


def test_populations():
    """
    >>> import microhapdb
    >>> p = microhapdb.populations
    >>> p[p.ID == 'SA000040E']
               ID     Name  Source
    48  SA000040E  Kachari  ALFRED
    >>> p[p.ID == 'SA000936S']
               ID     Name  Source
    53  SA000936S  Koreans  ALFRED
    >>> p[p.Name == 'Han']
               ID Name  Source
    29  SA004059S  Han  ALFRED
    30  SA004058R  Han  ALFRED
    31  SA000001B  Han  ALFRED
    32  SA000009J  Han  ALFRED
    >>> p.query('Name.str.contains("Afr")')
              ID               Name                        Source
    1     Africa             Africa  10.1016/j.fsigen.2018.05.008
    2  SA004047P  African Americans                        ALFRED
    3  SA000101C  African Americans                        ALFRED
    4  SA004242M    Afro-Caribbeans                        ALFRED
    """
    pop = microhapdb.populations
    assert pop.shape == (101, 3)
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
        ID Name  Source
 SA004059S  Han  ALFRED
 SA004058R  Han  ALFRED
 SA000001B  Han  ALFRED
 SA000009J  Han  ALFRED
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
'''
    terminal = capsys.readouterr()
    assert terminal.out.strip() == testout.strip()


@pytest.mark.parametrize('ident,data', [
    ('SA000019K', 'SA000019K  Russians  ALFRED'),
    ('Asia', 'Asia  Asia  10.1016/j.fsigen.2018.05.008'),
    ('Swedish', 'Swedish  Swedish  ISFG2019:P597'),
    ('HiroakiCohort', 'HiroakiCohort  HiroakiCohort  10.1016/j.legalmed.2015.06.003'),
])
def test_all_sources(ident, data, capsys):
    pop = microhapdb.populations[microhapdb.populations.ID == ident]
    microhapdb.population.print_table(pop)
    terminal = capsys.readouterr()
    print(terminal.out)
    assert data in terminal.out
