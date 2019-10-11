# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------


import builtins
from collections import defaultdict
from gzip import open as gzopen
import sys


def smartopen(filename, mode):
    """Smart file handler

    Determines whether to create a compressed or un-compressed file handle
    based on filename extension.
    """
    if mode not in ('r', 'w'):
        raise ValueError('invalid mode "{}"'.format(mode))
    if filename in ['-', None]:
        filehandle = sys.stdin if mode == 'r' else sys.stdout
        return filehandle
    openfunc = builtins.open
    if filename.endswith('.gz'):
        openfunc = gzopen
        mode += 't'
    return openfunc(filename, mode)


def parse_markers(instream):
    for line in instream:
        if not line.startswith('mh'):
            continue
        name, chrom, amplstart, amplend, localoffsets = line.strip().split(',')
        offsets = [int(amplstart) + int(o) for o in localoffsets.split(':')]
        offsets = sorted(offsets)
        if amplstart > amplend:
            amplstart, amplend = amplend, amplstart
        yield name, chrom, int(amplstart), int(amplend), offsets


def marker_variants(markerstream, vcfstream):
    rsids_by_pos = defaultdict(dict)
    for line in vcfstream:
        if line.startswith('#'):
            continue
        chrom, posstr, rsid, *values = line.split('\t')
        if not chrom.startswith('chr'):
            chrom = 'chr' + chrom
        pos = int(posstr) - 1
        rsids_by_pos[chrom][pos] = rsid

    for name, chrom, astart, aend, offsets in parse_markers(markerstream):
        rsids = list()
        for o in offsets:
            if o in rsids_by_pos[chrom]:
                rsid = rsids_by_pos[chrom][o]
                rsids.append(rsid)
            else:
                # Only variant not present in dbSNP
                assert (chrom, o) == ('chr22', 44857884), (chrom, o)
        yield name, chrom, offsets, rsids



#cli = argparse.ArgumentParser()
#cli.add_argument('lovd', type=argparse.FileType('r'))
#cli.add_argument('dbSNP', type=argparse.FileType('r'))
#args = cli.parse_args()
#
#class Locus(object):
#    def __init__(self, line):
#        v = line.strip().split(',')
#        self._values = v
#        self.label = v[0]
#        self.chrom = v[1]
#        self.start, self.end = int(v[2]), int(v[3])
#        self.positions = [int(o) + self.start for o in v[5].split(':')]
#        self.agg_alleles = dict()
#        self.alleles = [set() for _ in range(len(self.positions))]
#        self.rsids = defaultdict(list)
#
#    def add_allele(self, allelestr):
#        name, alleles = allelestr.strip().split(',')
#        self.agg_alleles[name] = [a.strip() for a in alleles.split(':')]
#        for n, a in enumerate(self.agg_alleles[name]):
#            self.alleles[n].add(a)
#
#    @property
#    def variants(self):
#        for position, alleleset in zip(self.positions, self.alleles):
#            yield self.chrom, position, alleleset
#
#    @property
#    def slug(self):
#        return '{:s}:{:d}-{:d}'.format(self.chrom, self.start, self.end)
#
#    def add_dbsnp_rsid(self, slug, rsid):
#        self.rsids[slug].append(rsid)
#
#
#loci = dict()
#locus = None
#for line in args.lovd:
#    if line.startswith('mh'):
#        locus = Locus(line)
#        loci[locus.label] = locus
#    else:
#        locus.add_allele(line)
#
#loci_by_position = dict()
#alleles_by_position = dict()
#for label, locus in loci.items():
#    for chrom, position, alleles in locus.variants:
#        slug = '{}:{}'.format(chrom, position)
#        alleles_by_position[slug] = alleles
#        loci_by_position[slug] = locus
#
#for line in args.dbSNP:
#    values = line.strip().split('\t')
#    chrom = 'chr' + values[0]
#    position = int(values[1]) - 1
#    slug = '{}:{}'.format(chrom, position)
#    rsID = values[2]
#    if slug not in alleles_by_position:
#        continue
#    locus = loci_by_position[slug]
#    lovd_alleles = alleles_by_position[slug]
#    refr = values[3].split(',')
#    alt = values[4].split(',')
#    dbsnp_alleles = set(refr + alt)
#    if lovd_alleles - dbsnp_alleles == set():
#        locus.add_dbsnp_rsid(slug, rsID)
#    elif slug == 'chr15:24802331':
#        locus.add_dbsnp_rsid(slug, 'rs12914037')  # doi:10.1016/j.fsigen.2018.05.008
#
#print('Label', 'Chrom', 'Position', 'Offset', 'Alleles', 'rsIDs', 'Source', sep='\t')
#for label, locus in loci.items():
#    for chrom, position, alleles in locus.variants:
#        slug = '{}:{}'.format(chrom, position)
#        source = 'dbSNP151'
#        rsids = '-'
#        if slug in locus.rsids:
#            rsids = ','.join(locus.rsids[slug])
#        else:
#            source = 'doi:10.1016/j.fsigen.2018.05.008'
#        alleles = ','.join(sorted(alleles))
#        print(locus.label, chrom, position, position - locus.start, alleles, rsids, source, sep='\t')
