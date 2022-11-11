# -------------------------------------------------------------------------------------------------
# Copyright (c) 2022, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from collections import Counter
from io import StringIO
import microhapdb
from warnings import warn


class Population:
    """Convenience class for accessing and manipulating population data

    >>> for pop in microhapdb.Population.from_ids(["CDX", "CHB", "CHS"]):
    ...   print(pop.popid, pop.name, pop.source)
    CDX Chinese Dai in Xishuangbanna, China 1KGP
    CHB Han Chinese in Beijing, China 1KGP
    CHS Southern Han Chinese 1KGP
    >>> for pop in microhapdb.Population.from_query("Name.str.contains('Japan')"):
    ...   print(pop)
    MHDBP-63967b883e	Japanese	10.1016/j.legalmed.2015.06.003
    SA000010B	Japanese	ALFRED
    JPT	Japanese in Tokyo, Japan	1KGP
    >>> pop = next(microhapdb.Population.from_ids(["SA004309Q"]))
    >>> print(pop.detail)
    --------------------------------------------------------------[ MicroHapDB ]----
    Iranian    (SA004309Q; source=ALFRED)
    - 399 total allele frequencies available
      for 65 markers
    # Alleles | # Markers
    ---------------------
            19|*
            14|*
            13|*
            12|****
             9|**
             8|*********
             7|***
             6|****
             5|**********
             4|******************************
    --------------------------------------------------------------------------------
    >>> microhapdb.populations.shape
    (109, 3)
    """
    def __init__(self, popid, name, source):
        self.popid = popid
        self.name = name
        self.source = source

    @classmethod
    def table_from_ids(cls, identifiers):
        ids = cls.standardize_ids(identifiers)
        table = microhapdb.populations[microhapdb.populations.ID.isin(ids)]
        return table

    @classmethod
    def table_from_query(cls, query):
        table = microhapdb.populations.query(query, engine="python")
        return table

    @classmethod
    def from_id(cls, identifier):
        table = cls.table_from_ids([identifier])
        if len(table) < 1:
            raise ValueError(f"popualtion {identifier} not found")
        if len(table) > 1:
            warn(f"ambiguous population identifier {identifier}", UserWarning)
        entry = table.iloc[0]
        return cls(entry.ID, entry.Name, entry.Source)

    @classmethod
    def from_ids(cls, identifiers):
        table = cls.table_from_ids(identifiers)
        yield from cls.objectify(table)

    @classmethod
    def from_query(cls, query):
        table = cls.table_from_query(query)
        yield from cls.objectify(table)

    @classmethod
    def objectify(cls, table):
        for i, row in table.iterrows():
            population = cls(row.ID, row.Name, row.Source)
            yield population

    @staticmethod
    def standardize_ids(identifiers):
        xref_ids = microhapdb.idmap[microhapdb.idmap.Xref.isin(identifiers)].ID
        ids = set(xref_ids) | set(identifiers)
        result = microhapdb.populations[
            (microhapdb.populations.ID.isin(ids)) | (microhapdb.populations.Name.isin(ids))
        ]
        return sorted(result.ID)

    def __str__(self):
        return "\t".join((self.popid, self.name, self.source))

    @property
    def detail(self):
        result = microhapdb.frequencies[microhapdb.frequencies.Population == self.popid]
        markers_with_n_alleles = Counter(result.groupby("Marker").size())
        output = StringIO()
        print(
            "--------------------------------------------------------------[ MicroHapDB ]----",
            file=output,
        )
        print(f"{self.name}    ({self.popid}; source={self.source})\n", file=output)
        print(f"- {len(result.Frequency)} total allele frequencies available", file=output)
        print(f"  for {len(result.Marker.unique())} markers\n", file=output)
        print("# Alleles | # Markers\n---------------------", file=output)
        for nalleles, nmarkers in sorted(markers_with_n_alleles.items(), reverse=True):
            print(f"       {nalleles:3d}|{'*' * nmarkers}", file=output)
        print(
            "--------------------------------------------------------------------------------",
            file=output,
        )
        return output.getvalue()
