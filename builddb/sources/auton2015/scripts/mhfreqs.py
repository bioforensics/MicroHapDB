#!/usr/bin/env python
# -------------------------------------------------------------------------------------------------
# Copyright (c) 2019, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from argparse import ArgumentParser
from collections import defaultdict
from glob import glob
import json
import os
import pandas
import rsidx
import sqlite3
import sys


def cli():
    """Command line interface"""
    parser = ArgumentParser()
    parser.add_argument(
        "--configfile", default="../../config.json", help="path to file containing database config"
    )
    parser.add_argument(
        "samplepops",
        help="two-column table with 1KGP sample in the first column and the corresponding 1KGP 3-character population code in the second column",
    )
    parser.add_argument(
        "markerrsids",
        help="two-column table with marker identifier in the first column and a comma-separated list of corresponding rsIDs in the second column",
    )
    return parser


def load_population_data(filename):
    """Load the population code for each 1KGP sample"""
    popsdf = pandas.read_csv(filename, sep="\t")
    return pandas.Series(popsdf.Population.values, index=popsdf.Sample).to_dict()


def construct_haplotypes(rsids, vcf, rsidxfile):
    """Construct haplotypes from each marker's list of rsIDS"""
    samples = list()
    haplo1 = defaultdict(list)
    haplo2 = defaultdict(list)
    with sqlite3.connect(rsidxfile) as dbconn:
        for line in rsidx.search.search(rsids, dbconn, vcf, header=True):
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    samples = line.strip().split()[9:]
                continue
            fields = line.strip().split()
            rsid = fields[2]
            if rsid not in rsids:
                continue
            ref = fields[3]
            alt = fields[4].split(",")
            samplegts = fields[9:]
            for sample, gt in zip(samples, samplegts):

                def _getallele(allelelabel):
                    if allelelabel == "0":
                        return ref
                    else:
                        i = int(allelelabel) - 1
                        return alt[i]

                if "|" in gt:
                    # Diploid autosomes
                    for altpresent, haplotype in zip(gt.split("|"), (haplo1, haplo2)):
                        allele = _getallele(altpresent)
                        haplotype[sample].append(allele)
                else:
                    # Haploid sex chromosome
                    allele = _getallele(gt)
                    haplo1[sample].append(allele)
    return samples, haplo1, haplo2


def compute_pop_counts(samples, haplo1, haplo2, samplepops):
    """Tally up observed haplotypes by population"""
    popallelecounts = defaultdict(lambda: defaultdict(int))
    for sample in samples:
        population = samplepops[sample]
        hap1gt = ",".join(haplo1[sample])
        popallelecounts[population][hap1gt] += 1
        if sample in haplo2:  # males don't have a second haplotype for chr0X
            hap2gt = ",".join(haplo2[sample])
            popallelecounts[population][hap2gt] += 1
    return popallelecounts


def mhfreqs(marker, rsids, vcffile, rsidxfile, samplepops):
    """Estimate microhaplotype frequences from the 1KGP data

    Construct haplotypes for each marker, compute a tally of allelic
    combinations observed in each population, and then compute the frequency
    for each.
    """
    samples, haplo1, haplo2 = construct_haplotypes(rsids, vcffile, rsidxfile)
    popallelecounts = compute_pop_counts(samples, haplo1, haplo2, samplepops)
    numvar_counts = defaultdict(int)
    exp_numvars = len(rsids)
    for population, allelecounts in popallelecounts.items():
        total = sum(allelecounts.values())
        for allele, count in sorted(allelecounts.items()):
            obs_numvars = allele.count(",") + 1
            numvar_counts[obs_numvars] += 1
            if obs_numvars != exp_numvars:
                continue
            freq = "{:.3f}".format(count / total)
            yield population, allele.replace(",", "|"), freq
    assert len(numvar_counts) == 1
    obs_numvars, count = list(numvar_counts.items())[0]
    if obs_numvars != exp_numvars:
        print(
            f"WARNING: unexpected number of variants found for marker {marker}:",
            f"expected {exp_numvars}, found {count} haplotypes with {obs_numvars} variants;",
            "WARNING: discarding haplotypes with unexpected number of variants",
            file=sys.stderr,
        )


def get_indexes_for_marker(marker, path):
    chrom = marker[2:4]
    if chrom[0] == "0":
        chrom = chrom[1]
    prefix = os.path.join(path, f"chr{chrom}")
    vcf = glob(f"{prefix}.vcf.gz")[0]
    rsidx = glob(f"{prefix}.rsidx")[0]
    return vcf, rsidx


def main(args):
    """Driver function"""
    with open(args.configfile, "r") as fh:
        config = json.load(fh)
    samplepops = load_population_data(args.samplepops)
    marker_rsids = dict()
    with open(args.markerrsids, "r") as fh:
        next(fh)
        for line in fh:
            marker, rsids = line.strip().split("\t")
            marker_rsids[marker] = rsids.split(",")
    print("Marker", "Population", "Allele", "Frequency", sep=",")
    for marker, rsids in marker_rsids.items():
        vcf, rsidx = get_indexes_for_marker(marker, path=config["1kgp_dir"])
        if not os.path.exists(vcf) or not os.path.exists(rsidx):
            raise RuntimeError("VCF and/or RSIDX files do not exist")
        for values in mhfreqs(marker, rsids, vcf, rsidx, samplepops):
            print(marker, *values, sep=",")


if __name__ == "__main__":
    main(cli().parse_args())
