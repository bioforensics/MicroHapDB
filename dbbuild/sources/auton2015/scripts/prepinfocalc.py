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
import json
from mhfreqs import cli, construct_haplotypes, get_indexes_for_marker, load_population_data


def main(args):
    with open(args.configfile, "r") as fh:
        config = json.load(fh)
    samplepops = load_population_data(args.samplepops)
    marker_rsids = dict()
    with open(args.markerrsids, "r") as fh:
        next(fh)
        for line in fh:
            marker, rsids = line.strip().split("\t")
            marker_rsids[marker] = rsids.split(",")

    # Print header: marker identifiers
    markers = sorted(marker_rsids.keys())
    print(*markers)

    # Load haplotypes into memory
    marker_haplotypes = defaultdict(dict)
    marker_haplotype_codes = defaultdict(dict)
    code = 0
    for marker in markers:
        vcf, rsidx = get_indexes_for_marker(marker, path=config["1kgp_dir"])
        rsids = marker_rsids[marker]
        samples, haplo1, haplo2 = construct_haplotypes(rsids, vcf, rsidx)

        def _get_hap_code(haplotype):
            nonlocal code
            if haplotype in marker_haplotype_codes[marker]:
                return marker_haplotype_codes[marker][haplotype]
            else:
                code += 1
                marker_haplotype_codes[marker][haplotype] = code
                return code

        for sample, gtlist in haplo1.items():
            haplotype = ",".join(gtlist)
            hapcode = _get_hap_code(haplotype)
            marker_haplotypes[sample][marker] = [hapcode]
        for sample, gtlist in haplo2.items():
            haplotype = ",".join(gtlist)
            hapcode = _get_hap_code(haplotype)
            marker_haplotypes[sample][marker].append(hapcode)

    # Print out each individual's genotype
    for sample in samples:
        data = [sample, samplepops[sample]]
        hap1codes = [marker_haplotypes[sample][mk][0] for mk in markers]
        print(*data, ".", ".", ".", *hap1codes)
        hap2codes = list()
        for mk in markers:
            code = -1
            if len(marker_haplotypes[sample][mk]) == 2:
                code = marker_haplotypes[sample][mk][1]
            hap2codes.append(code)
        print(*data, ".", ".", ".", *hap2codes)


if __name__ == "__main__":
    main(cli().parse_args())
