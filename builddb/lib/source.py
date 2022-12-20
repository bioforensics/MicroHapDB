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
import rsidx
import sqlite3
import subprocess
from tempfile import TemporaryDirectory


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


class VariantIndex:
    def __init__(self, table, dbsnp_path, chain_path):
        self.table = table
        self.dbsnp_path = Path(dbsnp_path)
        self.chain_path = Path(chain_path)
        self.coords_by_rsid = dict(GRCh37=dict(), GRCh38=dict())
        self.position_mapping = dict(GRCh37=dict(), GRCh38=dict())
        self.resolve_rsids()
        self.resolve_positions()

    @staticmethod
    def table_from_filenames(filenames):
        return pd.concat([pd.read_csv(fn) for fn in filenames]).reset_index()

    def resolve_rsids(self):
        rsids = set(self.all_rsids())
        vcf = self.dbsnp_path / "dbSNP_GRCh37.vcf.gz"
        idx = self.dbsnp_path / "dbSNP_GRCh37.rsidx"
        VariantIndex.rsidx_search(rsids, vcf, idx, self.coords_by_rsid["GRCh37"])
        vcf = self.dbsnp_path / "dbSNP_GRCh38.vcf.gz"
        idx = self.dbsnp_path / "dbSNP_GRCh38.rsidx"
        VariantIndex.rsidx_search(rsids, vcf, idx, self.coords_by_rsid["GRCh38"])

    def all_rsids(self):
        for n, row in self.table.iterrows():
            if pd.isna(row.VarRef):
                continue
            for rsid in row.VarRef.split(";"):
                yield rsid

    @staticmethod
    def rsidx_search(rsids, vcf, idx, vardict):
        with sqlite3.connect(idx) as db:
            for line in rsidx.search.search(rsids, db, vcf):
                values = line.strip().split("\t")
                pos = int(values[1])
                rsid = values[2]
                assert ";" not in rsid
                vardict[rsid] = pos

    def resolve_positions(self):
        self.map_positions("GRCh37", self.chain_path / "hg19ToHg38.over.chain.gz")
        self.map_positions("GRCh38", self.chain_path / "hg38ToHg19.over.chain.gz")

    def all_positions(self, refr):
        for n, row in self.table.iterrows():
            if pd.isna(row.Refr) or row.Refr != refr:
                continue
            for position in row.Positions.split(";"):
                yield row.Chrom, int(position) - 1, int(position)

    def map_positions(self, refr, chain):
        with TemporaryDirectory() as tmpdir:
            source_file = Path(tmpdir) / f"{refr}-source.bed"
            dest_file = Path(tmpdir) / f"{refr}-dist.bed"
            unmapped_file = Path(tmpdir) / f"{refr}-unmapped.bed"
            source_positions = pd.DataFrame(self.all_positions(refr), columns=["Chrom", "Start", "End"])
            source_positions.to_csv(source_file, sep="\t", header=False, index=False)
            args = map(str, ["liftOver", source_file, chain, dest_file, unmapped_file])
            subprocess.run(args)
            assert unmapped_file.stat().st_size == 0
            dest_positions = pd.read_csv(dest_file, sep="\t", names=["Chrom", "Start", "End"])
            for source, dest in zip(source_positions.iterrows(), dest_positions.iterrows()):
                self.position_mapping[refr][source[1].End] = dest[1].End














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
