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

from collections import defaultdict
import pandas as pd
from pyfaidx import Fasta as FastaIdx
from pysam import VariantFile


class GenotypeIndexer:
    def __init__(self, vcf, refr_fasta):
        self.vcf = VariantFile(vcf)
        self.refrseq_index = FastaIdx(refr_fasta)

    def get_haplotypes(self, chrom, positions):
        reference = self.get_reference_haplotype(chrom, positions)
        haplotypes = dict()
        for sample in self.samples:
            haplotypes[sample] = (Haplotype(reference), Haplotype(reference))
        for record in self.vcf.fetch(chrom, min(positions) - 1, max(positions)):
            if record.pos not in positions:
                continue
            for sample in self.samples:
                for haplotype, allele in zip(haplotypes[sample], record.samples[sample].alleles):
                    haplotype[record.pos] = allele
        for sample in self.samples:
            yield sample, haplotypes[sample][0], haplotypes[sample][1]

    @property
    def samples(self):
        return self.vcf.header.samples

    def get_reference_haplotype(self, chrom, positions):
        reference_haplotype = dict()
        for position in sorted(positions):
            allele = self.refrseq_index[chrom][position - 1]
            reference_haplotype[position] = allele
        return reference_haplotype


class Haplotype:
    def __init__(self, reference_positions):
        self.alleles = {pos: allele for pos, allele in reference_positions.items()}

    def __len__(self):
        return len(self.alleles)

    def __str__(self):
        return "|".join(map(str, self.alleles.values()))

    def __getitem__(self, position):
        if position not in self.alleles:
            raise IndexError(f"position {position} not set for this haplotype")
        return self.alleles[position]

    def __setitem__(self, position, allele):
        if position not in self.alleles:
            raise IndexError(f"position {position} not set for this haplotype")
        self.alleles[position] = allele
