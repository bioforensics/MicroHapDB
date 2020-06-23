# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from collections import defaultdict
from io import StringIO
from math import ceil
import microhapdb
from microhapdb.retrieve import id_in_series
import pandas


class TargetAmplicon():
    def __init__(self, marker, delta=10, minlen=80):
        if isinstance(marker, str):
            markerids = standardize_ids([marker])
            assert len(markerids) == 1
            marker = microhapdb.markers[microhapdb.markers.Name.isin(markerids)].iloc[0]
        self.data = marker
        self.delta = delta
        self.minlen = minlen
        self.flanking_sequence_data = microhapdb.sequences[
            microhapdb.sequences.Marker == self.data.Name
        ].iloc[0]

    @property
    def xrefs(self):
        xreflist = [self.data.PermID] + sorted(
            microhapdb.idmap[microhapdb.idmap.ID == self.data.Name].Xref
        )
        return xreflist

    @property
    def varrefs(self):
        return microhapdb.variantmap[microhapdb.variantmap.Marker == self.data.Name].Variant

    @property
    def offsets(self):
        return list(map(int, self.data.Offsets.split(',')))

    @property
    def variant_lengths(self):
        nvars = self.alleles[0].count(',') + 1
        lengths = [1] * nvars
        for n, row in microhapdb.indels[microhapdb.indels.Marker == self.data.Name].iterrows():
            assert row.VariantIndex < nvars, (row.VariantIndex, nvars)
            varalleles = [row.Refr] + row.Alt.split(',')
            varallelelengths = [len(va) for va in varalleles]
            lengths[row.VariantIndex] = max(varallelelengths)
        return lengths

    @property
    def reference_lengths(self):
        lengths = list()
        ind = microhapdb.indels
        for n in range(len(self.alleles[0].split(','))):
            refrlength = 1
            result = ind[(ind.Marker == self.data.Name) & (ind.VariantIndex == n)]
            if len(result) == 1:
                refrlength = len(result.Refr.iloc[0])
            lengths.append(refrlength)
        return lengths

    @property
    def slug(self):
        seqid = self.data.Chrom
        offsets = self.offsets
        return '{:s}:{:d}-{:d}'.format(seqid, offsets[0], offsets[-1])

    @property
    def marker_extent(self):
        o = self.offsets
        vl = self.variant_lengths
        return min(o), max(o) + vl[-1]

    @property
    def amplicon_interval(self):
        o = self.offsets
        vl = self.variant_lengths
        start, end = min(o) - self.delta, max(o) + self.delta + vl[-1]
        length = end - start
        if length < self.minlen:
            diff = self.minlen - length
            extension = ceil(diff / 2)
            start -= extension
            end += extension
        return start, end

    @property
    def start(self):
        return self.amplicon_interval[0]

    @property
    def end(self):
        return self.amplicon_interval[1]

    def __len__(self):
        return self.end - self.start

    @property
    def amplicon_offsets(self):
        start, end = self.amplicon_interval
        return [o - start for o in self.offsets]

    @property
    def marker_offsets(self):
        offsets = self.offsets
        return [o - offsets[0] for o in offsets]

    @property
    def alleles(self):
        return sorted(
            microhapdb.frequencies[microhapdb.frequencies.Marker == self.data.Name].Allele.unique()
        )

    @property
    def marker_seq(self):
        fullseq = self.flanking_sequence_data.Sequence
        mstart, mend = self.marker_extent
        seqstart = mstart - self.flanking_sequence_data.LeftFlank
        seqend = seqstart + (mend - mstart)
        markerseq = fullseq[seqstart:seqend]
        return markerseq

    @property
    def defline(self):
        varstring = ','.join(map(str, self.amplicon_offsets))
        line = '{name:s} PermID={pid:s} GRCh38:{slug:s} variants={var:s}'.format(
            name=self.data.Name, pid=self.data.PermID, slug=self.slug, var=varstring
        )
        result = microhapdb.idmap[microhapdb.idmap.ID == self.data.Name]
        if len(result) > 0:
            xrefstr = ','.join(result.Xref)
            line += ' Xref={x:s}'.format(x=xrefstr)
        return line

    @property
    def amplicon_seq(self):
        fullseq = self.flanking_sequence_data.Sequence
        astart, aend = self.amplicon_interval
        seqstart = astart - self.flanking_sequence_data.LeftFlank
        seqend = seqstart + (aend - astart)
        ampseq = fullseq[seqstart:seqend]
        return ampseq

    @property
    def fasta(self):
        out = StringIO()
        print('>', self.defline, sep='', file=out)
        ampseq = self.amplicon_seq
        if len(ampseq) < 80:
            print(ampseq, file=out)
        else:
            i = 0
            while i < len(ampseq):
                print(ampseq[i:i + 80], file=out)
                i += 80
        return out.getvalue().strip()

    def global_to_local(self, coord):
        if coord < self.start or coord > self.end:
            return None
        return coord - self.start

    def local_to_global(self, coord):
        if coord >= len(self):
            return None
        return coord + self.start

    def print_detail_definition(self, out):
        print('Marker Definition (GRCh38)', file=out)
        estart, eend = self.marker_extent
        extent = '{chr:s}:{start:d}-{end:d}'.format(chr=self.data.Chrom, start=estart, end=eend)
        extentstr = '    Marker extent\n        - {ext:s} ({length:d} bp)'.format(
            ext=extent, length=eend - estart
        )
        print(extentstr, file=out)
        astart, aend = self.amplicon_interval
        amplicon = '{chr:s}:{start:d}-{end:d} ({length:d} bp)'.format(
            chr=self.data.Chrom, start=astart, end=aend, length=aend - astart
        )
        print('    Target amplicon\n        -', amplicon, file=out)
        print('    Constituent variants', file=out)
        print('        - chromosome offsets:', self.data.Offsets, file=out)
        print(
            '        - marker offsets:',
            ','.join([str(o) for o in self.marker_offsets]),
            file=out
        )
        print(
            '        - amplicon offsets:',
            ','.join([str(o) for o in self.amplicon_offsets]),
            file=out
        )
        print('        - cross-references:', ', '.join(self.varrefs), file=out)
        print('    Observed haplotypes', file=out)
        for allele in self.alleles:
            print('        -', allele, file=out)
        print('\n', file=out)

    def print_detail_markerseq(self, out):
        print('--[ Core Marker Sequence ]--\n>', self.data.Name, sep='', file=out)
        markerseq = self.marker_seq
        if len(markerseq) < 80:
            print(markerseq, file=out)
        else:
            i = 0
            while i < len(markerseq):
                print(markerseq[i:i + 80], file=out)
                i += 80
        print('\n', file=out)

    def print_detail_ampliconseq(self, out):
        print('--[ Target Amplicon Sequence with Haplotypes ]--', file=out)
        blocks = list()
        prev = 0
        zipper = zip(self.amplicon_offsets, self.variant_lengths, self.reference_lengths)
        for n, (o, vl, rl) in enumerate(zipper):
            o_prev = o - prev
            if o_prev > 0:
                blocks.append(('span', o_prev))
            blocks.append(('variant', vl))
            prev = o + rl
        final = len(self.amplicon_seq) - o - vl
        if final > 0:
            blocks.append(('span', final))

        # Top row: variant indicators
        for blocktype, blocklength in blocks[:-1]:
            char = ' ' if blocktype == 'span' else '*'
            print(char * blocklength, end='', file=out)
        if blocks[-1][0] == 'variant':
            # only print final block if it's a variant; eliminate trailing whitespace
            print('*' * blocks[-1][1], end='', file=out)
        print('', file=out)

        # Second row: amplicon sequence
        i = 0
        n = -1
        refrlengths = self.reference_lengths
        for blocktype, blocklength in blocks:
            if blocktype == 'span':
                print(self.amplicon_seq[i:i + blocklength], end='', file=out)
                i += blocklength
            else:
                n += 1
                refrseq = '{seq:-<{length:d}s}'.format(
                    seq=self.amplicon_seq[i:i + refrlengths[n]],
                    length=blocklength,
                )
                print(refrseq, end='', file=out)
                i += refrlengths[n]
        print('', file=out)

        # Row 3+: alleles
        for allele in self.alleles:
            allelevars = allele.split(',')
            n = -1
            for blocktype, blocklength in blocks:
                if blocktype == 'span':
                    print('.' * blocklength, end='', file=out)
                else:
                    n += 1
                    allelestr = '{allele:-<{length:d}s}'.format(
                        allele=allelevars[n], length=self.variant_lengths[n]
                    )
                    print(allelestr, end='', file=out)
            print('', file=out)

    def __str__(self):
        out = StringIO()
        print('-' * 62, '[ MicroHapDB ]----', sep='', file=out)
        print(self.data.Name, '   a.k.a', ', '.join(self.xrefs), end='\n\n', file=out)
        self.print_detail_definition(out)
        self.print_detail_markerseq(out)
        self.print_detail_ampliconseq(out)
        print('-' * 80, file=out)
        return out.getvalue()


def standardize_ids(idents):
    ids = set()
    for ident in idents:
        if id_in_series(ident, microhapdb.variantmap.Variant):
            markernames = microhapdb.variantmap[microhapdb.variantmap.Variant == ident].Marker
            ids.update(markernames)
        elif id_in_series(ident, microhapdb.markers.PermID):
            result = microhapdb.markers[microhapdb.markers.PermID == ident]
            ids.update(result.Name)
        elif id_in_series(ident, microhapdb.markers.Name):
            result = microhapdb.markers[microhapdb.markers.Name == ident]
            ids.update(result.Name)
        elif id_in_series(ident, microhapdb.idmap.Xref):
            markername = microhapdb.idmap[microhapdb.idmap.Xref == ident].ID
            assert len(markername) == 1
            markername = markername.iloc[0]
            result = microhapdb.markers[microhapdb.markers.Name == markername]
            ids.update(result.Name)
    return microhapdb.markers[microhapdb.markers.Name.isin(ids)].Name


def print_table(table, delta=None, minlen=None, trunc=True):
    if trunc is not True:
        colwidth = pandas.get_option('display.max_colwidth')
        pandas.set_option('display.max_colwidth', 1000000)
    print(table.to_string(index=False))
    if trunc is not True:
        pandas.set_option('display.max_colwidth', colwidth)


def print_fasta(table, delta=10, minlen=80, trunc=None):
    for n, row in table.iterrows():
        amplicon = TargetAmplicon(row, delta=delta, minlen=minlen)
        print(amplicon.fasta)


def print_detail(table, delta=10, minlen=80, trunc=None):
    for n, row in table.iterrows():
        amplicon = TargetAmplicon(row, delta=delta, minlen=minlen)
        print(amplicon)
