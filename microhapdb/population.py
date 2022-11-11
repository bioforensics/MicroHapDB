# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from collections import Counter
from io import StringIO
import microhapdb


class Population:
    def __init__(self, popid, name, source):
        self.popid = popid
        self.name = name
        self.source = source

    @classmethod
    def from_ids(cls, identifiers):
        ids = cls.standardize_ids(identifiers)
        table = microhapdb.populations[microhapdb.populations.ID.isin(ids)]
        return table

    @classmethod
    def from_query(cls, query):
        table = microhapdb.populations.query(query, engine="python")
        return table

    @classmethod
    def objs_from_ids(cls, identifiers):
        table = cls.from_ids(identifiers)
        yield from cls.objectify(cls, table)

    @classmethod
    def objs_from_query(cls, query):
        table = cls.from_query(query)
        yield from cls.objectify(cls, table)

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
