# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from collections import defaultdict
from math import ceil
import microhapdb
from io import StringIO


class TargetAmplicon():
    def __init__(self, marker, delta=25, minlen=250):
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
    def marker_extent(self):
        o = self.offsets
        return min(o), max(o) + 1

    @property
    def amplicon_interval(self):
        o = self.offsets
        start, end = min(o) - self.delta, max(o) + self.delta + 1
        length = end - start
        if length < self.minlen:
            diff = self.minlen - length
            extension = ceil(diff / 2)
            start -= extension
            end += extension
        return start, end

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
    def amplicon_seq(self):
        fullseq = self.flanking_sequence_data.Sequence
        astart, aend = self.amplicon_interval
        seqstart = astart - self.flanking_sequence_data.LeftFlank
        seqend = seqstart + (aend - astart)
        ampseq = fullseq[seqstart:seqend]
        return ampseq

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
        print('    Observed alleles', file=out)
        for allele in self.alleles:
            print('        -', allele, file=out)
        print('\n', file=out)

    def print_detail_markerseq(self, out):
        print('--[ Marker Sequence ]--\n>', self.data.Name, sep='', file=out)
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
        print('--[ Target Amplicon Sequence with Alleles ]--', file=out)
        lengths = defaultdict(int)
        for allele in self.alleles:
            suballeles = allele.split(',')
            for n, a in enumerate(suballeles):
                if len(a) > lengths[n]:
                    lengths[n] = len(a)
        lengths = [lengths[l] for l in range(len(suballeles))]

        prev = 0
        for o, l in zip(self.amplicon_offsets, lengths):
            o_prev = o - prev
            print(' ' * o_prev, '*' * l, sep='', end='', file=out)
            prev = o + l
        print('', file=out)

        print(self.amplicon_seq, sep='', file=out)
        for allele in self.alleles:
            prev = 0
            for o, a, l in zip(self.amplicon_offsets, allele.split(','), lengths):
                o_prev = o - prev
                print(
                    '.' * o_prev, '{allele:-<{length:d}s}'.format(allele=a, length=l), sep='',
                    end='', file=out
                )
                prev = o + l
            final = len(self.amplicon_seq) - o - 1
            print('.' * final, file=out)

    def __str__(self):
        out = StringIO()
        print('-' * 58, '[ MicroHapulator ]---', file=out)
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


def print_table(table, **kwargs):
    print(table.to_string(index=False))


def print_detail(table, delta=25, minlen=250):
    for n, row in table.iterrows():
        amplicon = TargetAmplicon(row, delta=delta, minlen=minlen)
        print(amplicon)
