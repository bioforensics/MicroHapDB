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


class MicrohaplotypeContext:
    accessions = {
        "chr1": "NC_000001.11",
        "chr2": "NC_000002.12",
        "chr3": "NC_000003.12",
        "chr4": "NC_000004.12",
        "chr5": "NC_000005.10",
        "chr6": "NC_000006.12",
        "chr7": "NC_000007.14",
        "chr8": "NC_000008.11",
        "chr9": "NC_000009.12",
        "chr10": "NC_000010.11",
        "chr11": "NC_000011.10",
        "chr12": "NC_000012.12",
        "chr13": "NC_000013.11",
        "chr14": "NC_000014.9",
        "chr15": "NC_000015.10",
        "chr16": "NC_000016.10",
        "chr17": "NC_000017.11",
        "chr18": "NC_000018.10",
        "chr19": "NC_000019.10",
        "chr20": "NC_000020.11",
        "chr21": "NC_000021.9",
        "chr22": "NC_000022.11",
        "chrX": "NC_000023.11",
    }
    chromosomes = {accession: chromosome for chromosome, accession in accessions.items()}

    def __init__(self, marker):
        self.marker = marker

    def slug(self, delta=100):
        chrom, interval = self.marker.slug.split(":")
        accession = self.accessions[chrom]
        start, end = interval.split("-")
        if delta:
            start = int(start) - delta
            end = int(end) + delta
        return f"{accession}:{start}-{end}"
