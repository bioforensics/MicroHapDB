#!/usr/bin/env python

from argparse import ArgumentParser
from allel import GenotypeArray, weir_cockerham_fst
from collections import defaultdict
import numpy as np
import pandas
import sys


def load_sample_genotypes(instream):
    firstline = next(instream).strip()
    markers = firstline.split()
    pop_samples = defaultdict(set)
    genotypes = defaultdict(lambda: defaultdict(list))
    for line in instream:
        sample, population, _1, _2, _3, *alleles = line.strip().split()
        pop_samples[population].add(sample)
        for i, allele in enumerate(alleles):
            genotypes[sample][markers[i]].append(int(allele))
    return genotypes, pop_samples, markers


def build_genotype_array(genotypes, pop_samples, markers):
    g = list()
    for marker in markers:
        mgt = list()
        for population, sample_list in pop_samples.items():
            for sample in sorted(sample_list):
                mgt.append(genotypes[sample][marker])
        g.append(mgt)
    gt = GenotypeArray(g)
    print(
        'GenotypeArray construction complete:', gt.n_variants, 'markers,',
        gt.n_samples, 'samples, and a ploidy of', gt.ploidy, file=sys.stderr
    )
    return gt


def build_subpop_array(pop_samples):
    subpops = list()
    offset = 0
    for population, sample_list in pop_samples.items():
        subpops.append([offset, offset + len(sample_list) - 1])
        offset += len(sample_list)
    return subpops


def calculate_per_marker_fst(gt, subpops):
    a, b, c = weir_cockerham_fst(gt, subpops)
    fst = (
        np.sum(a, axis=1)
        /
        (np.sum(a, axis=1) + np.sum(b, axis=1) + np.sum(c, axis=1))
    )
    return fst


def main():
    cli = ArgumentParser()
    cli.add_argument('infile')
    args = cli.parse_args()

    with open(args.infile, 'r') as fh:
        genotypes, pop_samples, markers = load_sample_genotypes(fh)
    gt = build_genotype_array(genotypes, pop_samples, markers)
    subpops = build_subpop_array(pop_samples)
    fst = calculate_per_marker_fst(gt, subpops)
    data = pandas.DataFrame({
        'Marker': markers,
        'AvgFst': fst,
    })
    data.round({'AvgFst': 4}).to_csv(sys.stdout, sep='\t', index=False)


if __name__ == '__main__':
    main()
