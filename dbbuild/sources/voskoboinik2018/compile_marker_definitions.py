#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse
import os
import pandas
import re
import subprocess
import sys


class Variant():
    def __init__(self, line):
        self._fields = line.strip().split('\t')

    def __str__(self):
        return '\t'.join(self._fields)

    @property
    def rsid(self):
        return self._fields[2]

    @property
    def location(self):
        return int(self._fields[0]), int(self._fields[1])

    @property
    def af_slug(self):
        slug = list()
        for kvp in self._fields[7].split(';'):
            if kvp.count('=') != 1:
                continue
            key, value = kvp.split('=')
            if key.endswith('AF'):
                slug.append(key + '=' + value)
        return ';'.join(slug)

    @property
    def is_snp(self):
        for allele in self._fields[3].split(',') + self._fields[4].split(','):
            if len(allele) > 1:
                return False
        return True

    def is_rare(self, mincount=2):
        countmatch = re.search(r'AC=(\d+)', self._fields[7])
        assert countmatch
        count = int(countmatch.group(1))
        if count < mincount:
            return True
        return False


def parse_vcf(infile):
    for line in infile:
        if line.startswith('#'):
            continue
        yield Variant(line)


def parse_marker_definitions(table1file):
    data = pandas.read_csv(table1file, sep='\t')
    data['End'] = data['GRCh37Pos'] + data['Length']
    for n, row in data.iterrows():
        name = f'mh{row.Chrom:02d}LV-{n+1:02d}'
        region = f'{row.Chrom}:{row.GRCh37Pos}-{row.End}'
        yield name, region


def get_variants(markername, region, path):
    chrom = int(markername[2:4])
    infile = os.path.join(path, f'chr{chrom:d}.vcf.gz')
    outfile = '{mn:s}-unfiltered.vcf'.format(mn=markername)
    command = ['tabix', infile, region]
    with open(outfile, 'w') as fh:
        subprocess.check_call(command, universal_newlines=True, stdout=fh)


def filter_variants(markername):
    infile = '{mn:s}-unfiltered.vcf'.format(mn=markername)
    outfile = '{mn:s}-filtered.vcf'.format(mn=markername)
    variants = dict()
    with open(infile, 'r') as infh, open(outfile, 'w') as outfh:
        for var in parse_vcf(infh):
            if not var.is_snp:
                continue
            if var.is_rare():
                continue
            if var.af_slug in variants:
                print(
                    'DUPLICATE', variants[var.af_slug].rsid, var.rsid, var.af_slug,
                    file=sys.stderr
                )
            else:
                variants[var.af_slug] = var
        variants_to_keep = sorted(variants.values(), key=lambda v: v.location)
        rsidmap = {
            'rs113012024': 'rs10987426'
        }
        for var in variants_to_keep:
            if var.rsid in rsidmap:
                var._fields[2] = rsidmap[var.rsid]
            print(var, file=outfh)


def transfer_variants(markername, vcf, rsidx, refr='GRCh38'):
    infile = '{mn:s}-filtered.vcf'.format(mn=markername)
    rsids = list()
    with open(infile, 'r') as fh:
        for line in fh:
            if line.strip() == '' or line.startswith('#'):
                continue
            values = line.split('\t')
            rsid = values[2]
            rsids.append(rsid)
    outfile = f'{markername}-{refr}.vcf'
    command = ['rsidx', 'search', vcf, rsidx] + rsids
    with open(outfile, 'w') as fh:
        subprocess.check_call(command, universal_newlines=True, stdout=fh)


def vcf_to_rsid_offsets(vcfstream):
    offsets = list()
    rsids = list()
    for line in vcfstream:
        chrom, pos, rsid, *values = line.strip().split('\t')
        markerchrom = 'chr' + chrom
        offsets.append(int(pos) - 1)
        rsids.append(rsid)
    return markerchrom, rsids, offsets


def finalize_marker_definitions(markernames):
    data = {
        'Name': list(),
        'Xref': list(),
        'NumVars': list(),
        'Chrom': list(),
        'OffsetsHg37': list(),
        'OffsetsHg38': list(),
        'VarRef': list(),
    }
    for name in markernames:
        infile = '{mn:s}-GRCh38.vcf'.format(mn=name)
        offsets = list()
        rsids = list()
        with open(f'{name}-GRCh37.vcf', 'r') as fh:
            chrom, rsids37, offsets37 = vcf_to_rsid_offsets(fh)
        with open(f'{name}-GRCh38.vcf', 'r') as fh:
            chrom, rsids38, offsets38 = vcf_to_rsid_offsets(fh)
        assert sorted(rsids37) == sorted(rsids38), (sorted(rsids37), sorted(rsids38))
        data['Name'].append(name)
        data['Xref'].append(None)
        data['NumVars'].append(len(rsids38))
        data['Chrom'].append(chrom)
        data['OffsetsHg37'].append(','.join(map(str, offsets37)))
        data['OffsetsHg38'].append(','.join(map(str, offsets38)))
        data['VarRef'].append(','.join(rsids38))
    df = pandas.DataFrame(data)
    df.to_csv('marker.tsv', sep='\t', index=False)


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('dbsnp37', help='dbSNP VCF for GRCh37/hg19')
    parser.add_argument('dbsnp38', help='dbSNP VCF for GRCh38')
    parser.add_argument('rsidx37', help='rsidx index for GRCh37/hg19')
    parser.add_argument('rsidx38', help='rsidx index for GRCh38')
    parser.add_argument('dir1kgp', help='directory containing 1000 Genomes Project VCFs')
    parser.add_argument('table1', help='Table 1 from Voskoboinik et al 2018')
    return parser


def main(args):
    markernames = list()
    for name, region in parse_marker_definitions(args.table1):
        markernames.append(name)
        get_variants(name, region, args.dir1kgp)
        filter_variants(name)
        transfer_variants(name, args.dbsnp37, args.rsidx37, refr='GRCh37')
        transfer_variants(name, args.dbsnp38, args.rsidx38, refr='GRCh38')
    finalize_marker_definitions(markernames)


if __name__ == '__main__':
    main(cli().parse_args())
