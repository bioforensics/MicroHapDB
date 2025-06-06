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

from .interval import IntervalIndex
from .locus import Locus
from .marker import Marker
from .variant import VariantIndex
from collections import defaultdict
from io import StringIO
import json
import pandas as pd
from pathlib import Path


class DataSource:
    def __init__(self, source_path, variant_index):
        self.path = Path(source_path)
        self.metadata = DataSource.meta_from_json(self.path / "source.json")
        self.indels = DataSource.data_from_csv(self.path / "indels.csv")
        self.frequencies = DataSource.data_from_csv(self.path / "frequency.csv")
        if self.frequencies is not None:
            self.frequencies["Source"] = self.name
        self.populations = DataSource.data_from_csv(self.path / "population.csv")
        if self.populations is not None:
            self.populations.drop(columns=["Xref"], inplace=True)
            self.populations["Source"] = self.name
        self.markers = None
        markerpath = self.path / "marker.csv"
        if markerpath.is_file():
            self.markers = list(Marker.from_csv(markerpath, variant_index, source=self))

    @staticmethod
    def meta_from_json(metapath):
        if not metapath.is_file():
            raise FileNotFoundError(metapath)
        with open(metapath, "r") as fh:
            metadata = json.load(fh)
            assert "name" in metadata
            assert "year" in metadata
            assert "doi" in metadata
            assert "description" in metadata
            return metadata

    @staticmethod
    def data_from_csv(csvpath):
        if not csvpath.is_file():
            return None
        return pd.read_csv(csvpath)

    @property
    def name(self):
        return self.metadata["name"]

    @property
    def year(self):
        return self.metadata["year"]

    @property
    def order(self):
        if "order" in self.metadata:
            return self.metadata["order"]
        return 0

    @property
    def sortkey(self):
        return self.year, self.order, self.name.lower()

    @property
    def description(self):
        return self.metadata["description"]

    @property
    def num_snps(self):
        if self.markers is None:
            return 0
        num_vars = sum([m.numvars for m in self.markers])
        return num_vars - self.num_indels

    @property
    def num_indels(self):
        if self.indels is None:
            return 0
        return len(self.indels)

    @property
    def num_frequency_populations(self):
        if self.frequencies is None:
            return 0
        return len(self.frequencies.Population.unique())

    @property
    def num_distinct_haplotypes(self):
        if self.frequencies is None:
            return 0
        return len(self.frequencies.groupby(["Marker", "Allele"]))

    @property
    def num_frequencies(self):
        if self.frequencies is None:
            return 0
        return len(self.frequencies)

    def rename_markers(self, name_map):
        for oldname, newname in name_map.items():
            print(f"[{self.name}] {oldname} --> {newname}")
        if self.indels is not None:
            self.indels.replace(name_map, inplace=True)
        if self.frequencies is not None:
            self.frequencies.replace(name_map, inplace=True)

    def __str__(self):
        output = StringIO()
        print(f"[{self.name}]", file=output)
        if self.markers is not None:
            indelstr = f" and {self.num_indels} indels" if self.num_indels > 0 else ""
            varstr = f"based on {self.num_snps} SNPs{indelstr}"
            print(f"  - {len(self.markers)} marker definitions {varstr}", file=output)
        if self.populations is not None:
            samplestr = "samples" if len(self.populations) > 1 else "sample"
            print(f"  - {len(self.populations)} population {samplestr}", file=output)
        if self.frequencies is not None:
            if self.num_frequency_populations == 1:
                numhaps = self.num_distinct_haplotypes
                message = f"  - {numhaps} haplotype frequencies in 1 population"
            else:
                message = f"  - frequencies for {self.num_distinct_haplotypes} distinct haplotypes in {self.num_frequency_populations} populations; {self.num_frequencies} total frequencies"
            print(message, file=output)
        return output.getvalue().strip()


class SourceIndex:
    def __init__(self, rootdir, dbsnp_path, chain_path, exclude=[]):
        self.rootdir = Path(rootdir)
        self.dbsnp_path = Path(dbsnp_path)
        self.chain_path = Path(chain_path)
        self.exclusion_list = exclude
        self.populate_variants()
        self.populate_sources()
        self._markers = list()
        self.interval_index = IntervalIndex()

    def populate_variants(self):
        csvs = self.rootdir.glob("*/marker.csv")
        table = VariantIndex.table_from_filenames(csvs)
        self.variant_index = VariantIndex(table, self.dbsnp_path, self.chain_path)

    def populate_sources(self):
        self.sources = list()
        for sourcepath in self.rootdir.iterdir():
            if not sourcepath.is_dir():
                continue
            source = DataSource(sourcepath, self.variant_index)
            if source.name not in self.exclusion_list:
                self.sources.append(source)

    def all_markers(self):
        for source in self.sources:
            if source.markers is not None:
                for marker in source.markers:
                    yield marker

    def update_marker_names(self):
        markers_by_locus = defaultdict(Locus)
        for marker in self.all_markers():
            if marker.name in self.interval_index.mergeables:
                marker.name = self.interval_index.mergeables[marker.name]
            markers_by_locus[marker.locus].add(marker)
        source_name_map = defaultdict(dict)
        for name, locus in markers_by_locus.items():
            for marker in locus.resolve():
                self._markers.append(marker)
            for sourcename, namedict in locus.source_name_map.items():
                source_name_map[sourcename].update(namedict)
        for source in sorted(self.sources, key=lambda s: s.sortkey):
            source.rename_markers(source_name_map[source.name])

    def interval_check(self):
        for marker in self.all_markers():
            self.interval_index.add(marker)
        self.interval_index.check()

    @property
    def markers(self):
        table = list()
        for marker in sorted(self._markers, key=lambda m: m.sortkey):
            table.append(marker.fields)
        return pd.DataFrame(table, columns=Marker.field_names)

    @property
    def indels(self):
        table = pd.concat([source.indels for source in self.sources if source.indels is not None])
        table = table.sort_values(["Marker", "VariantIndex"]).reset_index(drop=True)
        return table

    @property
    def frequencies(self):
        table = pd.concat(
            [source.frequencies for source in self.sources if source.frequencies is not None]
        )
        table["Allele"] = table["Allele"].str.replace("|", ":")
        table = table.sort_values(["Marker", "Population"]).reset_index(drop=True)
        return table

    @property
    def populations(self):
        table = pd.concat([source.populations for source in self.sources if source.populations is not None])
        table = table.sort_values(["Name", "Source"]).reset_index(drop=True)
        return table

    @property
    def merges(self):
        table = pd.DataFrame(self.interval_index.mergeables.items(), columns=["Derivative", "Original"])
        table = table.sort_values(["Original", "Derivative"]).reset_index(drop=True)
        return table

    def __str__(self):
        output = StringIO()
        for source in sorted(self.sources, key=lambda s: s.sortkey):
            print(source, file=output)
        return output.getvalue()
