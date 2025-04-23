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

from mh_context import MicrohaplotypeContext
from collections import defaultdict
from intervaltree import IntervalTree
from subprocess import run
from tqdm import tqdm


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
            if (
                variant.is_indel
                and variant.has_frequency_data
                and max(variant.alternate_frequencies) >= min_freq
                and variant.rsid not in self.rsids
            ):
                self.index[variant.chrom][variant.start : variant.end] = variant

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
        return self.index[chrom][offset - distance : offset + distance]


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
            raise ValueError(
                f"no GnomAD, dbGaP, or 1000 Genomes frequency data: {self.rsid} {attributes}"
            )
        frequencies = freq_data.split(":")[1].split(",")
        frequencies = ["0.0" if f == "." else f for f in frequencies]
        frequencies = [float(freq) for freq in frequencies]
        assert len(frequencies) == 1 + len(self.alts)
        return sorted(frequencies, reverse=True)

    @property
    def alternate_frequencies(self):
        return self.frequency[1:]
