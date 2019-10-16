# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from math import ceil
import microhapdb


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


def marker_view(data, delta=25, minlen=250):
    print('-----------------------------------------------------------[ MicroHapulator ]---')
    xrefs = [data.PermID] + sorted(microhapdb.idmap[microhapdb.idmap.ID == data.Name].Xref)
    print(data.Name, '   a.k.a', ', '.join(xrefs), end='\n\n')
    print('Marker Definition (GRCh38)')
    offsets = list(map(int, data.Offsets.split(',')))
    estart, eend = min(offsets), max(offsets) + 1
    extent = '{chr:s}:{start:d}-{end:d}'.format(chr=data.Chrom, start=estart, end=eend)
    print('    Marker extent\n        - {ext:s} ({l:d} bp)'.format(ext=extent, l=eend - estart))
    astart, aend = min(offsets) - delta, max(offsets) + delta + 1
    alength = aend - astart
    if alength < minlen:
        diff = minlen - alength
        extend = ceil(diff / 2)
        astart -= extend
        aend += extend
        alength = aend - astart
    amplicon = '{chr:s}:{start:d}-{end:d} ({l:d} bp)'.format(
        chr=data.Chrom, start=astart, end=aend, l=aend - astart
    )
    print('    Target amplicon\n        -', amplicon)
    ampoffsets = [o - astart for o in offsets]
    markeroffsets = [o - offsets[0] for o in offsets]
    varrefs = microhapdb.variantmap[microhapdb.variantmap.Marker == data.Name].Variant
    print('    Constituent variants')
    print('        - chromosome offsets:', data.Offsets)
    print('        - marker offsets:', ','.join([str(o) for o in markeroffsets]))
    print('        - amplicon offsets:', ','.join([str(o) for o in ampoffsets]))
    print('        - cross-references:', ', '.join(varrefs))
    alleles = sorted(
        microhapdb.frequencies[microhapdb.frequencies.Marker == data.Name].Allele.unique()
    )
    print('    Observed alleles')
    for allele in alleles:
        print('        -', allele)
    print('\n')

    print('--[ Amplicon Sequence ]--')
    prev = 0
    for o in ampoffsets:
        o_prev = o - prev
        print(' ' * o_prev, '*', sep='', end='')
        prev = o + 1
    print('')
    sequencedata = microhapdb.sequences[microhapdb.sequences.Marker == data.Name].iloc[0]
    fullseq = sequencedata.Sequence
    seqstart = astart - sequencedata.LeftFlank
    seqend = seqstart + alength
    ampseq = fullseq[seqstart:seqend]
    print(ampseq, sep='')
    for allele in alleles:
        prev = 0
        for o, a in zip(ampoffsets, allele.split(',')):
            o_prev = o - prev
            print('.' * o_prev, a, sep='', end='')
            prev = o + 1
        final = len(ampseq) - o - 1
        print('.' * final)
    print('--------------------------------------------------------------------------------')


def print_detail(table, delta=25, minlen=250):
    for n, row in table.iterrows():
        marker_view(row, delta=delta, minlen=minlen)
