#!/usr/bin/env python

# -------------------------------------------------------------------------------------------------
# Copyright (c) 2024, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from dataclasses import dataclass
from itertools import combinations
import logging
import networkx as nx
from scipy.stats import pmean
from tqdm import tqdm
from util import marker_distance

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("[design.py]")


class LinkageGraph(list):
    @classmethod
    def populate(cls, markers, thresholds):
        markers = list(markers)
        chroms = sorted(set([mh.chrom for mh in markers]))
        graphs = cls()
        for chrom in chroms:
            cmarkers = [mh for mh in markers if mh.chrom == chrom]
            graph = LinkageGraph.populate_chromosome_graph(chrom, cmarkers, thresholds)
            if len(graph) > 0:
                graphs.append(graph)
        return graphs

    @staticmethod
    def populate_chromosome_graph(chrom, markers, thresholds):
        graph = nx.Graph()
        short = [mh for mh in markers if thresholds.is_short_candidate(mh)]
        long = [mh for mh in markers if thresholds.is_long_candidate(mh)]
        LinkageGraph.populate_top_ranked_markers(graph, short, thresholds.max_per_chrom_short)
        LinkageGraph.populate_top_ranked_markers(graph, long, thresholds.max_per_chrom_long)
        for mh1, mh2 in combinations(graph.nodes, 2):
            if mh1.chrom != mh2.chrom or marker_distance(mh1, mh2) > thresholds.ld_distance:
                graph.add_edge(mh1, mh2)
        return graph

    @staticmethod
    def populate_top_ranked_markers(graph, markers, limit):
        markers = sorted(markers, key=lambda mh: mh.data.Ae, reverse=True)
        markers = markers[:limit]
        if len(markers) > 0:
            graph.add_nodes_from(markers)

    @staticmethod
    def best_clique(cliques):
        best_ae = 0.0
        best_clique = None
        for clique in cliques:
            # In changing the sum to the power mean I was hoping to weight markers with larger Ae
            # scores more heavily. In practice it didn't seem to make much of a difference.
            # --DSS 2024-11-29
            # aggregate_ae = sum(mh.data.Ae for mh in clique)
            aggregate_ae = pmean([mh.data.Ae for mh in clique], 2)
            if aggregate_ae > best_ae:
                best_ae = aggregate_ae
                best_clique = clique
        return best_clique

    def design_panel(self):
        logger.info(f"Enumerating maximal cliques for {len(self)} graphs...")
        maximal_cliques_per_graph = [list(nx.find_cliques(graph)) for graph in self]
        logger.info(f"Determining optimal clique for {len(self)} graphs...")
        panel = list()
        for cliques in tqdm(maximal_cliques_per_graph):
            best_clique = self.best_clique(cliques)
            if best_clique is not None:
                panel.extend(best_clique)
        panel = sorted(panel, key=lambda mh: mh.name)
        return panel


@dataclass
class LinkageGraphThresholds:
    max_length_short: int = 99
    max_length_long: int = 259
    min_ae_short: float = 3.0
    min_ae_long: float = 5.0
    max_per_chrom_short: int = 8
    max_per_chrom_long: int = 8
    ld_distance: float = 10e6

    def is_short_candidate(self, microhap):
        return len(microhap) <= self.max_length_short and microhap.data.Ae >= self.min_ae_short

    def is_long_candidate(self, microhap):
        passes_length_criterion = self.max_length_short < len(microhap) <= self.max_length_long
        return passes_length_criterion and microhap.data.Ae > self.min_ae_long
