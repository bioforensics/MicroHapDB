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

from .marker import BaseMarkerDefinition, CompleteMarkerDefinition
from .resolver import Resolver
from collections import Counter, defaultdict
from io import StringIO
import json
import pandas as pd
from pathlib import Path


class DataSource:
    def __init__(self, source_path, dbsnp_path=None):
        self.path = Path(source_path)
        self.metadata = DataSource.meta_from_json(self.path / "source.json")
        self.indels = DataSource.data_from_csv(self.path / "indels.csv")
        self.frequencies = DataSource.data_from_csv(self.path / "frequency.csv")
        self.populations = DataSource.data_from_csv(self.path / "population.csv")
        self.markers = None
        markerpath = self.path / "marker.csv"
        if markerpath.is_file():
            self.markers = list(BaseMarkerDefinition.from_csv(markerpath, source=self))

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
    def __init__(self, rootdir, dbsnp_dir):
        self.sources = SourceIndex.populate(rootdir)
        self.all_rsids = set()
        for source in self.sources:
            if source.markers is None:
                continue
            for marker in source.markers:
                if marker.varref is not None:
                    self.all_rsids.update(marker.varref)
        self.resolver = Resolver(dbsnp_dir)
        self.resolver.resolve_rsids(self.all_rsids)
        self.markers = list()
        self.update_marker_names()

    @staticmethod
    def populate(rootdir):
        sources = list()
        for sourcepath in Path(rootdir).iterdir():
            if not sourcepath.is_dir():
                continue
            sources.append(DataSource(sourcepath))
        return sources

    def update_marker_names(self):
        markers_by_name = defaultdict(list)
        for source in self.sources:
            if source.markers is not None:
                for marker in source.markers:
                    marker = CompleteMarkerDefinition(marker, self.resolver)
                    if marker.is_match is False:
                        print(f"[{marker.base.source.name}::{marker.name}] WARNING: inferred positions do not match published positions")
                    markers_by_name[marker.name].append(marker)
        for name, markers in markers_by_name.items():
            name_by_positions = dict()
            distinct_definitions = set([m.posstr for m in markers])
            for marker in sorted(markers, key=lambda m: (m.base.source.year, m.base.name.lower())):
                if len(markers) > 1:
                    if marker.posstr in name_by_positions:
                        print(f"Marker {marker.base.name} as defined in {marker.base.source.name} is redundant")
                        continue
                    elif len(distinct_definitions) > 1:
                        newname = f"{marker.base.name}.v{len(name_by_positions) + 1}"
                        name_by_positions[marker.posstr] = newname
                    marker.base.name = newname
                self.markers.append(marker)

    def marker_definitions(self):
        table = list()
        for marker in sorted(self.markers, key=lambda m: (m.chrom_num, m.span, m.base.name)):
            entry = (marker.base.name, marker.base.chrom, marker.posstr, marker.base.source.name)
            table.append(entry)
        return pd.DataFrame(table, columns=["Name", "Chrom", "Positions", "Source"])

    def __str__(self):
        output = StringIO()
        for source in sorted(self.sources, key=lambda s: (s.year, s.name.lower())):
            print(source, file=output)
        return output.getvalue()
