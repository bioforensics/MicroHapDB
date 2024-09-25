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

from io import StringIO
from math import ceil
import microhapdb
import pandas as pd


class Marker:
    """Convenience class for accessing and manipulating marker data

    >>> marker = microhapdb.Marker.from_id("mh13KK-218.v1")
    >>> marker.slug
    'chr13:53486692-53486837'
    >>> len(marker)
    146
    >>> marker.offsets
    [53486691, 53486745, 53486756, 53486836]
    >>> marker.start
    53486691
    >>> marker.marker_offsets
    [0, 54, 65, 145]
    >>> marker.target_offsets
    [10, 64, 75, 155]
    >>> marker.definition
              Marker  Offset  Chrom  ChromOffset
    0  mh13KK-218.v1      10  chr13     53486691
    1  mh13KK-218.v1      64  chr13     53486745
    2  mh13KK-218.v1      75  chr13     53486756
    3  mh13KK-218.v1     155  chr13     53486836
    >>> print(marker.fasta)
    >mh13KK-218.v1 GRCh38:chr13:53486681-53486847 variants=10,64,75,155
    ATAGCACATTTCCAAGTTGTTCTAGTGAATTACTGAACTGGATAGGATTGTGGAAACCTGTGAATAATAGCTAGGTAGTC
    AGAAGACATGGTGCGCTGGGGATCCTCAAAGTGTGGCTGTTAACTGAAATGAAGGTACTCTTGTGGAGGACTGAGCCCTT
    AACATG
    """

    def __init__(self, marker, delta=10, minlen=80, extendmode=0):
        self.extendmode = int(extendmode)
        self.data = marker
        self.delta = delta
        self.minlen = minlen

    @staticmethod
    def table_from_ids(identifiers):
        ids = Marker.standardize_ids(identifiers)
        table = microhapdb.markers[microhapdb.markers.Name.isin(ids)]
        return table

    @staticmethod
    def table_from_query(query):
        table = microhapdb.markers.query(query, engine="python")
        return table

    @staticmethod
    def definitions_from_ids(identifiers, **kwargs):
        return pd.concat([marker.definition for marker in Marker.from_ids(identifiers, **kwargs)])

    @staticmethod
    def table_from_region(regionstr):
        markers = microhapdb.markers.copy()
        chrom, start, end = Marker.parse_regionstr(regionstr)
        query = f'Chrom == "{chrom}"'
        if start is not None:
            query += f" and (Start < {end}"
            query += f" and End > {start})"
        result = markers.query(query)
        return result

    @classmethod
    def from_id(cls, identifier, **kwargs):
        markerids = cls.standardize_ids([identifier])
        if len(markerids) < 1:
            raise ValueError(f"no such marker '{identifier}'")
        result = microhapdb.markers[microhapdb.markers.Name.isin(markerids)]
        assert len(result) == 1, identifier
        marker = result.iloc[0]
        return cls(marker, **kwargs)

    @classmethod
    def from_ids(cls, identifiers, **kwargs):
        table = cls.table_from_ids(identifiers)
        yield from cls.objectify(table, **kwargs)

    @classmethod
    def from_query(cls, query, **kwargs):
        table = cls.table_from_query(query)
        yield from cls.objectify(table, **kwargs)

    @classmethod
    def from_region(cls, region, **kwargs):
        table = cls.table_from_region(region)
        yield from cls.objectify(table, **kwargs)

    @classmethod
    def objectify(cls, table, **kwargs):
        for i, row in table.iterrows():
            marker = cls(row, **kwargs)
            yield marker

    @staticmethod
    def parse_regionstr(regionstr):
        """Retrieve chromosome name and coordinates from a region string

        Region string is expected to be in one of the two following formats: 'chr3'
        or 'chr3:1000000-5000000'.

        >>> Marker.parse_regionstr("chr12")
        ('chr12', None, None)
        >>> Marker.parse_regionstr("chr12:345-678")
        ('chr12', 345, 678)
        """
        chrom, start, end = None, None, None
        if ":" in regionstr:
            chrom, rng = regionstr.split(":")
            if rng.count("-") != 1:
                raise ValueError(f'cannot parse region "{regionstr}"')
            startstr, endstr = rng.split("-")
            start, end = int(startstr), int(endstr)
        else:
            chrom = regionstr
        return chrom, start, end

    @staticmethod
    def standardize_ids(idents):
        def id_in_series(ident, series):
            return ident in series.values

        ids = set()
        for ident in idents:
            locusnames = microhapdb.markers.Name.apply(lambda x: x.split(".")[0])
            if id_in_series(ident, microhapdb.variantmap.Variant):
                result = microhapdb.variantmap[microhapdb.variantmap.Variant == ident]
                ids.update(result.Marker)
            elif id_in_series(ident, microhapdb.markers.Name) or id_in_series(ident, locusnames):
                result = microhapdb.markers[microhapdb.markers.Name.str.contains(ident)]
                ids.update(result.Name)
            elif id_in_series(ident, microhapdb.merged.Derivative):
                result = microhapdb.merged[microhapdb.merged.Derivative.str.contains(ident)]
                assert len(result) == 1
                original = result.Original.iloc[0]
                result = microhapdb.markers[microhapdb.markers.Name.str.contains(original)]
                ids.update(result.Name)
        return sorted(microhapdb.markers[microhapdb.markers.Name.isin(ids)].Name)

    def __str__(self):
        return f"{self.name} ({self.slug})"

    @property
    def chrom(self):
        return self.data.Chrom

    @property
    def chrom_num(self):
        if self.chrom == "chrX":
            return 23
        return int(self.chrom[3:])

    @property
    def name(self):
        return self.data.Name

    @property
    def locus(self):
        return self.name.split(".")[0]

    @property
    def slug(self):
        seqid = self.chrom
        return f"{seqid}:{self.start+1}-{self.end}"

    def __len__(self):
        return self.end - self.start

    @property
    def nvar(self):
        return len(self.offsets)

    @property
    def start(self):
        return int(self.data.Start - 1)

    @property
    def end(self):
        return int(self.data.End - 1 + self.variant_lengths[-1])

    @property
    def target_slug(self):
        seqid = self.chrom
        start, end = self.target_interval
        return f"{seqid}:{start}-{end}"

    @property
    def target_length(self):
        start, end = self.target_interval
        return end - start

    @property
    def target_interval(self):
        o = self.offsets
        vl = self.variant_lengths
        start, end = min(o) - self.delta, max(o) + self.delta + vl[-1]
        length = end - start
        if length < self.minlen:
            diff = self.minlen - length
            if self.extendmode < 0:
                start -= diff
            elif self.extendmode > 0:
                end += diff
            else:
                extension = ceil(diff / 2)
                start -= extension
                end += extension
        return start, end

    @property
    def variant_lengths(self):
        nvars = len(self.offsets)
        lengths = [1] * nvars
        for n, row in microhapdb.indels[microhapdb.indels.Marker == self.name].iterrows():
            assert row.VariantIndex < nvars, (row.VariantIndex, nvars)
            varalleles = [row.Refr] + row.Alt.split(",")
            varallelelengths = [len(va) for va in varalleles]
            lengths[row.VariantIndex] = max(varallelelengths)
        return lengths

    @property
    def detail(self):
        output = StringIO()
        print("-" * 62, "[ MicroHapDB ]----", sep="", file=output)
        self.print_detail_all_names(output)
        self.print_detail_definition(output)
        self.print_detail_markerseq(output)
        self.print_detail_targetseq(output)
        print("-" * 80, file=output)
        return output.getvalue()

    @property
    def varrefs(self):
        if pd.isna(self.data.RSIDs):
            return []
        return self.data.RSIDs.split(";")

    def print_detail_all_names(self, out):
        names = microhapdb.merged[microhapdb.merged.Original == self.locus].Derivative.unique()
        print("MH Locus:", self.locus, end="", file=out)
        if len(names) > 0:
            names = ", ".join(sorted(names))
            print(f" (a.k.a. {names})", end="", file=out)
        print(f"\nMarker:   {self.name} (Source: {self.data.Source})\n", file=out)

    def print_detail_definition(self, out):
        marker_slug = f"{self.slug} ({len(self)} bp)"
        target_slug = f"{self.target_slug} ({self.target_length} bp)"
        offsets37 = ", ".join([str(o) for o in self.offsets37])
        offsets38 = ", ".join([str(o) for o in self.offsets])
        marker_offsets = ", ".join([str(o) for o in self.marker_offsets])
        target_offsets = ", ".join([str(o) for o in self.target_offsets])
        print(f"Marker Definition", file=out)
        print(f"    Marker extent\n        - {marker_slug}", file=out)
        print(f"    Target locus\n        - {target_slug}", file=out)
        print(f"    Constituent variants", file=out)
        print(f"        - chromosome offsets (GRCh37):", offsets37, file=out)
        print(f"        - chromosome offsets (GRCh38):", offsets38, file=out)
        print(f"        - marker offsets:", marker_offsets, file=out)
        print(f"        - target offsets:", target_offsets, file=out)
        print(f"        - cross-references:", ", ".join(self.varrefs), file=out)
        print(f"    Observed haplotypes", file=out)
        for allele in self.alleles:
            print("        -", allele, file=out)
        print("\n", file=out)

    def print_detail_markerseq(self, out):
        print("--[ Core Marker Sequence ]--\n>", self.name, sep="", file=out)
        markerseq = self.marker_seq
        if len(markerseq) < 80:
            print(markerseq, file=out)
        else:
            i = 0
            while i < len(markerseq):
                print(markerseq[i : i + 80], file=out)
                i += 80
        print("\n", file=out)

    def print_detail_targetseq(self, out):
        print("--[ Marker Target Sequence with MH alleles (haplotypes) ]--", file=out)
        blocks = self.build_target_seq_blocks()
        self.print_detail_targetseq_variants(blocks, out)  # Top row: variant indicators
        self.print_detail_targetseq_sequence(blocks, out)  # Second row: target sequence
        self.print_detail_targetseq_alleles(blocks, out)  # Row 3+: alleles

    def build_target_seq_blocks(self):
        blocks = list()
        previous = 0
        iterator = zip(self.target_offsets, self.variant_lengths, self.reference_lengths)
        for offset, varlen, reflen in iterator:
            offset_prev = offset - previous
            if offset_prev > 0:
                blocks.append(("span", offset_prev))
            blocks.append(("variant", varlen))
            previous = offset + reflen
        final = len(self.target_seq) - offset - varlen
        if final > 0:
            blocks.append(("span", final))
        return blocks

    def print_detail_targetseq_variants(self, blocks, out):
        for blocktype, blocklength in blocks[:-1]:
            char = " " if blocktype == "span" else "*"
            print(char * blocklength, end="", file=out)
        if blocks[-1][0] == "variant":
            # only print final block if it's a variant; eliminate trailing whitespace
            print("*" * blocks[-1][1], end="", file=out)
        print("", file=out)

    def print_detail_targetseq_sequence(self, blocks, out):
        i = 0
        n = -1
        refrlengths = self.reference_lengths
        for blocktype, blocklength in blocks:
            if blocktype == "span":
                print(self.target_seq[i : i + blocklength], end="", file=out)
                i += blocklength
            else:
                n += 1
                refrseq = "{seq:-<{length:d}s}".format(
                    seq=self.target_seq[i : i + refrlengths[n]],
                    length=blocklength,
                )
                print(refrseq, end="", file=out)
                i += refrlengths[n]
        print("", file=out)

    def print_detail_targetseq_alleles(self, blocks, out):
        for allele in self.alleles:
            allelevars = allele.split("|")
            if len(allelevars) != len(self.variant_lengths):
                message = f"num var mismatch ({self.name}): {len(allelevars)} vs {len(self.variant_lengths)}"
                raise ValueError(message)
            n = -1
            for blocktype, blocklength in blocks:
                if blocktype == "span":
                    print("." * blocklength, end="", file=out)
                else:
                    n += 1
                    allelestr = "{allele:-<{length:d}s}".format(
                        allele=allelevars[n], length=self.variant_lengths[n]
                    )
                    print(allelestr, end="", file=out)
            print("", file=out)

    @property
    def marker_offsets(self):
        offsets = self.offsets
        return [o - offsets[0] for o in offsets]

    @property
    def target_offsets(self):
        start, end = self.target_interval
        return [o - start for o in self.offsets]

    @property
    def offsets(self):
        return sorted([int(p) - 1 for p in self.data.Positions.split(";")])

    @property
    def offsets37(self):
        return sorted([int(p) - 1 for p in self.data.Positions37.split(";")])

    @property
    def marker_seq(self):
        start = self.start
        end = self.end
        return str(microhapdb.hg38[self.chrom][start:end])

    @property
    def target_seq(self):
        start, end = self.target_interval
        return str(microhapdb.hg38[self.chrom][start:end])

    @property
    def flank_seqs(self):
        markerseq = self.target_seq
        first_snp_idx = self.target_offsets[0]
        last_snp_idx = self.target_offsets[-1]
        left = markerseq[:first_snp_idx]
        right = markerseq[last_snp_idx:]
        return left, right

    @property
    def reference_lengths(self):
        lengths = list()
        ind = microhapdb.indels
        for n in range(self.nvar):
            refrlength = 1
            result = ind[(ind.Marker == self.name) & (ind.VariantIndex == n)]
            if len(result) == 1:
                refrlength = len(result.Refr.iloc[0])
            lengths.append(refrlength)
        return lengths

    @property
    def alleles(self):
        return sorted(
            microhapdb.frequencies[microhapdb.frequencies.Marker == self.name].Allele.unique()
        )

    @property
    def defline(self):
        varstring = ",".join(map(str, self.target_offsets))
        parts = [
            self.name,
            f"GRCh38:{self.target_slug}",
            f"variants={varstring}",
        ]
        line = " ".join(parts)
        return line

    @property
    def fasta(self):
        out = StringIO()
        print(">", self.defline, sep="", file=out)
        targetseq = self.target_seq
        if len(targetseq) < 80:
            print(targetseq, file=out)
        else:
            i = 0
            while i < len(targetseq):
                print(targetseq[i : i + 80], file=out)
                i += 80
        return out.getvalue().strip()

    @property
    def definition(self):
        variants = list()
        for offset, refr_offset in zip(self.target_offsets, self.offsets):
            variants.append((self.name, offset, self.chrom, refr_offset))
        return pd.DataFrame(variants, columns=["Marker", "Offset", "Chrom", "ChromOffset"])

    def global_to_local(self, coord):
        start, end = self.target_interval
        if coord < start or coord > end:
            return None
        return coord - start

    def local_to_global(self, coord):
        start, end = self.target_interval
        if coord >= end - start:
            return None
        return coord + start


class Locus:
    def __init__(self, markers=None):
        self.markers = list() if markers is None else markers

    @property
    def definition(self):
        variants = list()
        for marker in sorted(self.markers, key=lambda m: m.name):
            start, end = self.target_interval
            offsets = [o - start for o in marker.offsets]
            for offset, refr_offset in zip(offsets, marker.offsets):
                variants.append((marker.name, offset, marker.chrom, refr_offset))
        return pd.DataFrame(variants, columns=["Marker", "Offset", "Chrom", "ChromOffset"])

    @property
    def fasta(self):
        out = StringIO()
        print(">", self.defline, sep="", file=out)
        targetseq = self.target_seq
        if len(targetseq) < 80:
            print(targetseq, file=out)
        else:
            i = 0
            while i < len(targetseq):
                print(targetseq[i : i + 80], file=out)
                i += 80
        return out.getvalue().strip()

    @property
    def defline(self):
        parts = [self.name, f"GRCh38:{self.target_slug}"]
        for marker in sorted(self.markers, key=lambda m: m.name):
            start, end = self.target_interval
            offsets = [o - start for o in marker.offsets]
            varstring = ",".join(map(str, offsets))
            parts.append(f"{marker.name}={varstring}")
        line = " ".join(parts)
        return line

    @property
    def target_seq(self):
        start, end = self.target_interval
        return str(microhapdb.hg38[self.chrom][start:end])

    @property
    def name(self):
        if len(self.markers) == 0:
            return None
        return self.markers[0].locus

    @property
    def target_slug(self):
        seqid = self.chrom
        start, end = self.target_interval
        return f"{seqid}:{start}-{end}"

    @property
    def target_interval(self):
        start = float("Inf")
        end = 0
        for marker in self.markers:
            mstart, mend = marker.target_interval
            if mstart < start:
                start = mstart
            if mend > end:
                end = mend
        return start, end

    @property
    def chrom(self):
        if len(self.markers) == 0:
            return None
        return self.markers[0].chrom
