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

from collections import Counter, defaultdict
import pandas as pd
from pathlib import Path
from pyfaidx import Fasta as FastaIdx
from pysam import VariantFile


def collect_haplotypes(markers, sample_pops, refr_file, vcfdir):
    haplo_table = list()
    for n, marker in markers.iterrows():
        vcfs = list(Path(vcfdir).glob(f"*{marker.Chrom}.*.vcf.gz"))
        if len(vcfs) != 1:
            message = f"found {len(vcfs)} VCFs for {marker.Chrom} ($VCFDIR/*{marker.Chrom}.*.vcf.gz), expected 1"
            raise FileNotFoundError(message)
        gti = GenotypeIndexer(vcfs[0], refr_file)
        positions = list(map(int, marker.Positions.split(";")))
        for sample, hapl1, hapl2 in gti.get_haplotypes(marker.Chrom, positions):
            if sample not in sample_pops:
                continue
            pop, superpop, gender = sample_pops[sample]
            if marker.Chrom == "chrX" and gender == "male":
                hapl2 = None
            entry = (marker.Name, sample, pop, superpop, hapl1, hapl2)
            haplo_table.append(entry)
    colnames = ["Marker", "Sample", "Population", "Superpopulation", "Haplotype1", "Haplotype2"]
    haplotypes = pd.DataFrame(haplo_table, columns=colnames)
    return haplotypes


class GenotypeIndexer:
    def __init__(self, vcf, refr_fasta):
        self.vcf = VariantFile(vcf)
        self.refrseq_index = FastaIdx(refr_fasta)

    def get_haplotypes(self, chrom, positions):
        reference = self.get_reference_haplotype(chrom, positions)
        if len(reference) != len(positions):
            message = f"variant count mismatch: {len(reference)} vs {len(positions)}"
            raise ValueError(message)
        haplotypes = dict()
        for sample in self.samples:
            haplotypes[sample] = (Haplotype(reference), Haplotype(reference))
        for record in self.vcf.fetch(chrom, min(positions) - 1, max(positions)):
            if record.pos not in positions:
                continue
            for allele in record.alleles:
                if len(allele) > 1:
                    break
            else:
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


def compile_sample_populations(vcf_path, pop_table, pedigree, dofilters=True):
    vcf = VariantFile(vcf_path)
    samples = set(vcf.header.samples)
    superpops = dict(zip(pop_table["Population code"], pop_table["Superpopulation code"]))
    if dofilters:
        pedigree = pedigree[pedigree["Individual ID"].isin(samples)]
        pedigree = pedigree[~pedigree["Paternal ID"].isin(samples)]
        pedigree = pedigree[~pedigree["Maternal ID"].isin(samples)]
    pedigree["Superpopulation"] = pedigree["Population"].map(superpops)
    pedigree["Gender"] = pedigree["Gender"].map({1: "male", 2: "female"})
    pedigree = pedigree[["Individual ID", "Gender", "Population", "Superpopulation"]].rename(
        columns={"Individual ID": "Sample"}
    )
    summarize_population_groups(pedigree, superpops)
    return pedigree


def summarize_population_groups(pedigree, superpops):
    popcounts = Counter(pedigree["Population"])
    superpopcounts = defaultdict(dict)
    for pop, count in popcounts.items():
        superpop = superpops[pop]
        superpopcounts[superpop][pop] = count
    for superpop, counts in sorted(superpopcounts.items()):
        total = sum(counts.values())
        print(f"- {superpop}: {total}")
        for pop, count in sorted(counts.items()):
            print(f"  - {pop}: {count}")


def compile_frequencies(haplotypes):
    table = list()
    for freqdata in list_frequencies(haplotypes):
        table.append(freqdata)
    return pd.DataFrame(table, columns=["Marker", "Population", "Allele", "Frequency", "Count"])


def list_frequencies(haplotypes):
    admixed = ("ACB", "ASW", "CLM", "MXL", "PEL", "PUR")
    pop_tallies = defaultdict(lambda: defaultdict(Counter))
    for n, row in haplotypes.iterrows():
        for haplokey in ("Haplotype1", "Haplotype2"):
            mhallele = [row[haplokey]]
            if not pd.isna(mhallele):
                pop_tallies[row["Marker"]][row["Population"]].update(mhallele)
                if row["Population"] not in admixed:
                    pop_tallies[row["Marker"]][row["Superpopulation"]].update(mhallele)
    for marker, popcounts in sorted(pop_tallies.items()):
        for population, haplocounts in sorted(popcounts.items()):
            total_count = sum(haplocounts.values())
            for mhallele, count in sorted(haplocounts.items()):
                freq = count / total_count
                yield marker, population, mhallele, freq, total_count


def compute_aes(frequencies):
    aes = list()
    superpops = ("AFR", "AMR", "EAS", "EUR", "SAS")
    for marker, marker_data in frequencies.groupby("Marker"):
        population_aes = list()
        for population, pop_data in marker_data.groupby("Population"):
            ae = 1.0 / sum([f**2 for f in pop_data.Frequency])
            entry = (marker, population, ae)
            aes.append(entry)
            if population not in superpops:
                population_aes.append(ae)
        avg_ae = sum(population_aes) / len(population_aes)
        entry = (marker, "1KGP", avg_ae)
        aes.append(entry)
    return pd.DataFrame(aes, columns=["Marker", "Population", "Ae"])
