#!/usr/bin/env python

from argparse import ArgumentParser
from collections import defaultdict
import rsidx
import sqlite3


def cli():
    parser = ArgumentParser()
    parser.add_argument('vcf')
    parser.add_argument('rsidx')
    parser.add_argument('samplepops')
    parser.add_argument('rsids', nargs='+')
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
            alt = fields[4]
            samplegts = fields[9:]
            for sample, gt in zip(samples, samplegts):
                for altpresent, haplotype in zip(gt.split('|'), (haplo1, haplo2)):
                    allele = alt if altpresent == '1' else ref
                    haplotype[sample].append(allele)
    return samples, haplo1, haplo2


def compute_pop_counts(samples, haplo1, haplo2, samplepops):
    popallelecounts = defaultdict(lambda: defaultdict(int))
    for sample in samples:
        population = samplepops[sample]
        hap1gt = ','.join(haplo1[sample])
        hap2gt = ','.join(haplo2[sample])
        popallelecounts[population][hap1gt] += 1
        popallelecounts[population][hap2gt] += 1
    return popallelecounts


def mhfreqs(rsids, vcffile, rsidxfile, samplepopsfile):
    samplepops = load_population_data(samplepopsfile)
    samples, haplo1, haplo2 = construct_haplotypes(rsids, vcffile, rsidxfile)
    popallelecounts = compute_pop_counts(samples, haplo1, haplo2, samplepops)
    for population, allelecounts in popallelecounts.items():
        total = sum(allelecounts.values())
        for allele, count in allelecounts.items():
            freq = '{:.3f}'.format(count / total)
            yield population, allele, freq


def main(args):
    print('Population', 'Allele', 'Frequency', sep='\t')
    for values in mhfreqs(args.rsids, args.vcf, args.rsidx, args.samplepops):
        print(*values, sep='\t')


if __name__ == '__main__':
    main(cli().parse_args())
