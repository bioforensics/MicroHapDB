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

from itertools import combinations
import logging
import networkx as nx
from tqdm import tqdm

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("[design.py]")


class LinkageGraph(list):
    @classmethod
    def populate(
        cls, markers, batches=8, max_per_chrom_long=8, max_per_chrom_short=8, distance=10e6
    ):
        graphs = cls()
        while len(graphs) < batches:
            graphs.append(nx.Graph())
        markers = list(markers)
        chroms = sorted(set([mh.chrom for mh in markers]))
        for chrom in chroms:
            cmarkers = [mh for mh in markers if mh.chrom == chrom]
            short = [mh for mh in cmarkers if len(mh) < 100 and mh.data.Ae > 3.0]
            long = [mh for mh in cmarkers if 100 <= len(mh) < 260 and mh.data.Ae > 5.0]
            for marker_set, max_per_chrom in zip(
                (short, long), (max_per_chrom_short, max_per_chrom_long)
            ):
                marker_set = sorted(marker_set, key=lambda mh: mh.data.Ae, reverse=True)
                marker_set = marker_set[:max_per_chrom]
                if len(marker_set) == 0:
                    continue
                index = marker_set[0].chrom_num % batches
                graph = graphs[index]
                graph.add_nodes_from(marker_set)
        for n, graph in enumerate(graphs):
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
        logger.info(f"Enumerating maximal cliques for {len(self)} batches...")
        maximal_cliques_per_batch = [list(nx.find_cliques(graph)) for graph in tqdm(self)]
        logger.info(f"Determining optimal clique for {len(self)} batches...")
        panel = list()
        for n, cliques in enumerate(maximal_cliques_per_batch):
            logger.info(f"Batch {n+1}: comparing maximal cliques")
            best_clique = self.best_clique(cliques)
            panel.extend(best_clique)
        panel = sorted(panel, key=lambda mh: mh.name)
        return panel
