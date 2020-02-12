#!/usr/bin/env python

from argparse import ArgumentParser
from collections import defaultdict
from glob import glob
import os
import rsidx
import sqlite3


def cli():
    parser = ArgumentParser()
    parser.add_argument('samplepops')
    parser.add_argument('markerrsids')
    return parser


def load_population_data(filename):
    samplepops = dict()
    with open(filename, 'r') as fh:
        next(fh)
        for line in fh:
            sample, population = line.strip().split('\t')
            samplepops[sample] = population
    return samplepops


def construct_haplotypes(rsids, vcf, rsidxfile):
    samples = list()
    haplo1 = defaultdict(list)
    haplo2 = defaultdict(list)
    with sqlite3.connect(rsidxfile) as dbconn:
        for line in rsidx.search.search(rsids, dbconn, vcf, header=True):
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    samples = line.strip().split()[9:]
                continue
            fields = line.strip().split()
            ref = fields[3]
            alt = fields[4].split(',')
            samplegts = fields[9:]
            for sample, gt in zip(samples, samplegts):
                def _getallele(allelelabel):
                    if allelelabel == '0':
                        return ref
                    else:
                        i = int(allelelabel) - 1
                        return alt[i]
                if '|' in gt:
                    # Diploid autosomes
                    for altpresent, haplotype in zip(gt.split('|'), (haplo1, haplo2)):
                        allele = _getallele(altpresent)
                        haplotype[sample].append(allele)
                else:
                    # Haploid sex chromosome
                    allele = _getallele(gt)
                    haplo1[sample].append(allele)
    return samples, haplo1, haplo2


def compute_pop_counts(samples, haplo1, haplo2, samplepops):
    popallelecounts = defaultdict(lambda: defaultdict(int))
    for sample in samples:
        population = samplepops[sample]
        hap1gt = ','.join(haplo1[sample])
        popallelecounts[population][hap1gt] += 1
        if sample in haplo2:  # males don't have a second haplotype for chr0X
            hap2gt = ','.join(haplo2[sample])
            popallelecounts[population][hap2gt] += 1
    return popallelecounts


def mhfreqs(rsids, vcffile, rsidxfile, samplepops):
    samples, haplo1, haplo2 = construct_haplotypes(rsids, vcffile, rsidxfile)
    popallelecounts = compute_pop_counts(samples, haplo1, haplo2, samplepops)
    for population, allelecounts in popallelecounts.items():
        total = sum(allelecounts.values())
        for allele, count in sorted(allelecounts.items()):
            freq = '{:.3f}'.format(count / total)
            yield population, allele, freq


def main(args):
    samplepops = load_population_data(args.samplepops)

    marker_rsids = dict()
    with open(args.markerrsids, 'r') as fh:
        next(fh)
        for line in fh:
            marker, rsids = line.strip().split('\t')
            marker_rsids[marker] = rsids.split(',')

    print('Marker', 'Population', 'Allele', 'Frequency', sep='\t')
    for marker, rsids in marker_rsids.items():
        chrom = marker[2:4]
        if chrom[0] == '0':
            chrom = chrom[1]
        prefix = 'ALL.chr' + chrom + '.phase3_shapeit2_mvncall_integrated_v??.20130502.genotypes'
        vcf = glob(prefix + '.vcf.gz')[0]
        rsidx = glob(prefix + '.rsidx')[0]
        if not os.path.exists(vcf) or not os.path.exists(rsidx):
            raise RuntimeError('VCF and/or RSIDX files do not exist')
        for values in mhfreqs(rsids, vcf, rsidx, samplepops):
            print(marker, *values, sep='\t')


if __name__ == '__main__':
    main(cli().parse_args())
