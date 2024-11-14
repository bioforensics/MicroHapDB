#!/usr/bin/env python

from argparse import ArgumentParser
from collections import defaultdict
from itertools import combinations
import logging
import networkx as nx
import microhapdb
import polars as pl
import sys
from tqdm import tqdm

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("[design.py]")


def main2(
    rmsk_path, markers_to_keep=None, filtered_path=None, max_per_chrom=8, batches=8, distance=10e6
):
    markers = microhapdb.Marker.objectify(microhapdb.markers)
    markers = [mh for mh in markers if "LV" not in mh.name]
    if markers_to_keep is None:
        print("Parsing UCSC RepeatMasker annotations...", end="", file=sys.stderr)
        repeats = parse_ucsc_rmsk_track(rmsk_path)
        print("done!", file=sys.stderr)
        print("Filtering markers...", file=sys.stderr)
        markers = [mh for mh in filter_markers(markers, repeats)]
        if filtered_path:
            with open(filtered_path, "w") as fh:
                identifiers = sorted([mh.name for mh in markers])
                print(*identifiers, sep="\n", file=fh)
    else:
        with open(markers_to_keep, "r") as fh:
            keepers = set(fh.read().strip().split())
        markers = [mh for mh in markers if mh.name in keepers]
    print("Populating linkage graph...", file=sys.stderr)
    ###lgraph = LinkageGraph.populate(markers)
    lgraph = LinkageGraph2.populate(
        markers, max_per_chrom=max_per_chrom, batches=batches, distance=distance
    )
    # for n, graph in enumerate(lgraph):
    #    print(f"{n} Nodes: {graph.number_of_nodes()}, edges: {graph.number_of_edges()}", file=sys.stderr)
    cliques = list()
    for graph in tqdm(lgraph):
        cliques.append(list(nx.find_cliques(graph)))
    for n, clique_list in enumerate(cliques):
        print(f"{n} Cliques: {len(clique_list)}", file=sys.stderr)
    panel = lgraph.panel
    for mh in panel:
        print(mh.chrom, mh.name, len(mh), mh.data.Ae, file=sys.stderr)
    print(len(panel), sep="\n", file=sys.stderr)
    raise SystemExit(0)
    # for chrom, cliques in lgraph.cliques_by_chrom.items():
    #     counts = Counter([len(c) for c in cliques])
    #     print(chrom, *counts.most_common())
    panel = list()
    for chrom, clique in lgraph.best_cliques:
        panel.extend(clique)
    panel = sorted(panel, key=lambda mh: mh.name)
    for mh in panel:
        print(mh.chrom, mh.name, len(mh), mh.data.Ae)
    print(len(panel))


def main(masks, read_masked=None, write_masked=None, max_per_chrom=8, batches=8, distance=10e6):
    markers = microhapdb.Marker.objectify(microhapdb.markers)
    markers = list(markers)
    masking = MaskManager(markers, masks, premasked_path=read_masked)
    if write_masked:
        masking.write(write_masked)
    markers = [mh for mh in markers if mh.name in masking.identifiers]
    logger.info("Populating linkage graph")
    graph = LinkageGraph.populate(
        markers, max_per_chrom=max_per_chrom, batches=batches, distance=distance
    )
    panel = graph.panel
    for mh in panel:
        print(mh.chrom, mh.name, len(mh), mh.data.Ae, file=sys.stderr)
    print(len(panel), sep="\n", file=sys.stderr)


class LinkageGraphOriginal(defaultdict):
    @classmethod
    def populate(cls, markers):
        graph = cls(nx.Graph)
        markers = list(markers)
        chroms = sorted(set([mh.chrom for mh in markers]))
        for chrom in (pbar := tqdm(chroms)):
            pbar.set_description(chrom)
            cmarkers = [mh for mh in markers if mh.chrom == chrom]
            short = [mh for mh in cmarkers if len(mh) < 100 and mh.data.Ae > 3.0]
            long = [mh for mh in cmarkers if 100 <= len(mh) < 250 and mh.data.Ae > 5.0]
            for marker_set in (short, long):
                marker_set = sorted(marker_set, key=lambda mh: mh.data.Ae, reverse=True)
                marker_set = marker_set[:12]
                graph[chrom].add_nodes_from(marker_set)
            for mh1, mh2 in combinations(graph[chrom].nodes, 2):
                if abs(mh1.start - mh2.start) > 10e6:
                    graph[chrom].add_edge(mh1, mh2)
        return graph

    @property
    def cliques_by_chrom(self):
        cliques = dict()
        for chrom, graph in sorted(self.items(), key=lambda g: g[1].number_of_edges()):
            chrom_cliques = list(nx.find_cliques(graph, nodes=None))
            cliques[chrom] = chrom_cliques
        return cliques

    @property
    def best_cliques(self):
        chrom_cliques = self.cliques_by_chrom
        for chrom, cliques in chrom_cliques.items():
            best_ae = 0.0
            best_clique = None
            for clique in cliques:
                ae = sum(mh.data.Ae for mh in clique)
                if ae > best_ae:
                    best_ae = ae
                    best_clique = clique
            yield chrom, best_clique


class LinkageGraph(list):
    @classmethod
    def populate(cls, markers, batches=8, max_per_chrom=8, distance=10e6):
        graphs = cls()
        while len(graphs) < batches:
            graphs.append(nx.Graph())
        markers = list(markers)
        chroms = sorted(set([mh.chrom for mh in markers]))
        for chrom in (pbar := tqdm(chroms)):
            pbar.set_description(chrom)
            cmarkers = [mh for mh in markers if mh.chrom == chrom]
            short = [mh for mh in cmarkers if len(mh) < 100 and mh.data.Ae > 3.0]
            long = [mh for mh in cmarkers if 100 <= len(mh) < 260 and mh.data.Ae > 5.0]
            for marker_set in (short, long):
                marker_set = sorted(marker_set, key=lambda mh: mh.data.Ae, reverse=True)
                marker_set = marker_set[:max_per_chrom]
                if len(marker_set) == 0:
                    continue
                index = marker_set[0].chrom_num % batches
                graph = graphs[index]
                graph.add_nodes_from(marker_set)
        for graph in graphs:
            for mh1, mh2 in combinations(graph.nodes, 2):
                if mh1.chrom != mh2.chrom or abs(mh1.start - mh2.start) > distance:
                    graph.add_edge(mh1, mh2)
        return graphs

    @staticmethod
    def best_clique(cliques):
        best_ae = 0.0
        best_clique = None
        for clique in tqdm(cliques):
            ae = sum(mh.data.Ae for mh in clique)
            if ae > best_ae:
                best_ae = ae
                best_clique = clique
        return best_clique

    @property
    def panel(self):
        print(f"Enumerating maximal cliques for {len(self)} batches...", end="", file=sys.stderr)
        maximal_cliques_per_batch = [list(nx.find_cliques(graph)) for graph in tqdm(self)]
        print(f"Determining optimal clique for {len(self)} batches...", end="", file=sys.stderr)
        panel = list()
        for cliques in maximal_cliques_per_batch:
            best_clique = self.best_clique(cliques)
            panel.extend(best_clique)
        panel = sorted(panel, key=lambda mh: mh.name)
        return panel


class MaskManager:
    def __init__(self, markers, mask_paths, premasked_path=None):
        if premasked_path:
            logger.info("Loading pre-filtered markers")
            with open(premasked_path, "r") as fh:
                self.identifiers = set(fh.read().strip().split())
        else:
            logger.info(f"Filtering markers using {len(mask_paths)} user-provided mask(s)")
            masker = Masker(*mask_paths)
            self.markers = [mh for mh in masker.mask(markers)]
            loci = set([mh.locus for mh in self.markers])
            message = f"{len(self.identifiers)} ({len(loci)} loci) markers passed masking filter"
            logger.info(message)

    def write(self, path):
        logger.info(f"Writing filtered markers to {path}")
        with open(path, "w") as fh:
            print(*self.identifiers, sep="\n", file=fh)

    @property
    def identifiers(self):
        return sorted([mh.name for mh in self.markers])


class Masker:
    def __init__(self, *paths):
        columns = ["Chrom", "Start", "End"]
        tables = [
            pl.read_csv(path, sep="\t", new_columns=columns, has_header=False) for path in paths
        ]
        self.coords = pl.concat(tables)

    def mask(self, markers):
        for mh in (pbar := tqdm(markers)):
            pbar.set_description(f"{mh.name:<20}")
            overlap = self.coords.filter(
                (pl.col("Chrom") == mh.chrom)
                & (pl.col("End") > mh.start)
                & (pl.col("Start") < mh.end)
            )
            if len(overlap) == 0:
                yield mh


def get_parser():
    parser = ArgumentParser()
    masking = parser.add_argument_group("Masking")
    masking.add_argument(
        "--masks",
        nargs="+",
        metavar="PATH",
        help="path to one or more BED files containing genomic coordinates; microhaps overlapping these coordinates will be excluded from panel design",
    )
    masking.add_argument(
        "--write-masked",
        metavar="PATH",
        help="write identifiers of filtered markers to the given PATH",
    )
    masking.add_argument(
        "--read-masked",
        metavar="PATH",
        help="read marker identifiers from the given PATH of pre-filtered markers; disables masking if specified",
    )
    parser.add_argument(
        "--max-per-chrom",
        type=int,
        default=8,
        metavar="M",
        help="select the M highest ranked short and long markers per chromosome by Ae; by default M=8",
    )
    parser.add_argument(
        "--batches",
        type=int,
        default=8,
        metavar="B",
        help="split the linkage graph into B batches; by default B=8",
    )
    parser.add_argument(
        "--distance",
        type=float,
        default=10e6,
        metavar="D",
        help="two markers must be separated by more than D bp to be considered independently inherited (as a heuristic); by default D=10000000",
    )
    parser._optionals.title = "Miscellaneous"
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    main(
        args.masks,
        read_masked=args.read_masked,
        write_masked=args.write_masked,
        max_per_chrom=args.max_per_chrom,
        batches=args.batches,
        distance=args.distance,
    )
