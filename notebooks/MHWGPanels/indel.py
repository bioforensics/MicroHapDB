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

from argparse import ArgumentParser
from collections import defaultdict
from intervaltree import IntervalTree
import microhapdb
import pandas as pd
from subprocess import run
from tqdm import tqdm


def main(markers, dbsnp_path, delta=100, min_freq=0.005):
    index = IndelIndex()
    index.populate(markers, dbsnp_path, delta=delta, min_freq=min_freq)
    print(f"Found {len(index)} indels with an alt allele at frequency ≥ {min_freq}")
    for mh in markers:
        for offset in mh.offsets:
            indels = index.query(mh.chrom, offset, distance=4)
            if len(indels) > 0:
                indels = [interval.data for interval in indels]
                freq_data = "|".join([f"{indel.rsid}:{indel.frequency}" for indel in indels])
                max_freq = max([max(indel.alternate_frequencies) for indel in indels])
                print(mh.name, len(mh), mh.data.Ae, max_freq, offset, freq_data, sep="\t")
                break


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


class IndelIndex:
    def __init__(self):
        self.index = defaultdict(IntervalTree)
        self.rsids = set()

    def __len__(self):
        return sum([len(tree) for tree in self.index.values()])

    def populate(self, markers, dbsnp_path, delta=100, min_freq=0.005):
        for buffer in self.context_stream(markers):
            self.process_buffer(buffer, dbsnp_path, delta=delta, min_freq=min_freq)

    def process_buffer(self, buffer, dbsnp_path, delta=100, min_freq=0.005):
        regions = [context.slug(delta=delta) for context in buffer]
        process = run(["tabix", dbsnp_path, *regions], capture_output=True, text=True)
        for line in process.stdout.strip().split("\n"):
            variant = Variant(line)
            if variant.is_indel and variant.has_frequency_data and max(variant.alternate_frequencies) >= min_freq and variant.rsid not in self.rsids:
                self.index[variant.chrom][variant.start:variant.end] = variant

    @staticmethod
    def context_stream(marker_stream, size=25):
        buffer = list()
        for mh in (pb := tqdm(marker_stream)):
            pb.set_description(f"{mh.name:<20}")
            context = MicrohaplotypeContext(mh)
            buffer.append(context)
            if len(buffer) == size:
                yield buffer
                buffer = list()
        if len(buffer) > 0:
            yield buffer

    def query(self, chrom, offset, distance=4):
        return self.index[chrom][offset]


class Variant:
    def __init__(self, data):
        self.fields = data.split("\t")

    @property
    def chrom(self):
        return MicrohaplotypeContext.chromosomes[self.fields[0]]

    @property
    def is_indel(self):
        for alt in self.alts:
            if len(self.refr) != len(alt):
                return True
        else:
            return False

    @property
    def rsid(self):
        return self.fields[2]

    @property
    def refr(self):
        allele = self.fields[3]
        assert "," not in allele
        return allele

    @property
    def alts(self):
        return self.fields[4].split(",")

    @property
    def start(self):
        return int(self.fields[1]) - 1

    @property
    def end(self):
        allele_lengths = [len(allele) for allele in [[self.refr] + self.alts]]
        return int(self.fields[1]) + max(allele_lengths)

    @property
    def is_common(self):
        return "COMMON" in self.fields[7]

    @property
    def has_frequency_data(self):
        if "FREQ=" not in self.fields[7]:
            return False
        for source in ("GnomAD:", "dbGaP_PopFreq:", "1000Genomes:"):
            if source in self.fields[7]:
                return True
        else:
            return False

    @property
    def frequency(self):
        attributes = self.fields[7]
        freq = [attr for attr in attributes.split(";") if attr.startswith("FREQ=")]
        if len(freq) == 0:
            return None
        assert len(freq) == 1, (self.rsid, attributes)
        data = freq[0].split("=")[1]
        for source in ("GnomAD:", "dbGaP_PopFreq:", "1000Genomes:"):
            if source in data:
                freq_data = [entry for entry in data.split("|") if entry.startswith(source)]
                assert len(freq_data) == 1, (self.rsid, attributes)
                freq_data = freq_data[0]
                break
        else:
            raise ValueError(f"no GnomAD, dbGaP, or 1000 Genomes frequency data: {self.rsid} {attributes}")
        frequencies = freq_data.split(":")[1].split(",")
        frequencies = ["0.0" if f == "." else f for f in frequencies]
        frequencies = [float(freq) for freq in frequencies]
        assert len(frequencies) == 1 + len(self.alts)
        return frequencies

    @property
    def alternate_frequencies(self):
        freqs = sorted(self.frequency, reverse=True)
        return freqs[1:]


def get_parser():
    parser = ArgumentParser()
    parser.add_argument("markers")
    parser.add_argument("dbsnp_path")
    parser.add_argument("--delta", type=int, default=25, metavar="D")
    parser.add_argument("--aes")
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    marker_table = pd.read_csv(args.markers)
    marker_table["Ae"] = 0.0
    if args.aes:
        aes = pd.read_csv(args.aes)
        popaes = aes[aes.Population == "1KGP"].drop(columns=["Population"])
        marker_table = marker_table.drop(columns=["Ae"]).join(popaes.set_index("Marker"), on="Name")
    markers = list(microhapdb.Marker.objectify(marker_table))
    main(markers, args.dbsnp_path)
