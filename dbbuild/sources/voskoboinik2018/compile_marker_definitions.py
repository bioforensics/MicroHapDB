#!/usr/bin/env python
import argparse
import json
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


def get_regions(table1file):
    data = pandas.read_csv(table1file, sep='\t')
    data['End'] = data['GRCh37Pos'] + data['Length']
    for n, row in data.iterrows():
        yield '{chr:d}:{start:d}-{end:d}'.format(chr=row.Chrom, start=row.GRCh37Pos, end=row.End)


def get_marker_names(table1file):
    data = pandas.read_csv(table1file, sep='\t')
    for n, row in data.iterrows():
        yield 'mh{chr:02d}LV-{sn:02d}'.format(chr=row.Chrom, sn=n+1)


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


def transfer_variants(markername, vcf, rsidx):
    infile = '{mn:s}-filtered.vcf'.format(mn=markername)
    rsids = list()
    with open(infile, 'r') as fh:
        for line in fh:
            if line.strip() == '' or line.startswith('#'):
                continue
            values = line.split('\t')
            rsid = values[2]
            rsids.append(rsid)
    outfile = '{mn:s}-GRCh38.vcf'.format(mn=markername)
    command = ['rsidx', 'search', vcf, rsidx] + rsids
    with open(outfile, 'w') as fh:
        subprocess.check_call(command, universal_newlines=True, stdout=fh)


def finalize_marker_definitions(markernames):
    data = {
        'Name': list(),
        'Xref': list(),
        'Reference': list(),
        'Chrom': list(),
        'Offsets': list(),
        'VarRef': list(),
    }
    for name in markernames:
        infile = '{mn:s}-GRCh38.vcf'.format(mn=name)
        offsets = list()
        rsids = list()
        with open(infile, 'r') as fh:
            for line in fh:
                chrom, pos, rsid, *values = line.strip().split('\t')
                markerchrom = 'chr' + chrom
                offsets.append(int(pos) - 1)
                rsids.append(rsid)
        offsetstr = ','.join(map(str, offsets))
        rsidstr = ','.join(rsids)
        data['Name'].append(name)
        data['Xref'].append(None)
        data['Reference'].append('GRCh38')
        data['Chrom'].append(markerchrom)
        data['Offsets'].append(offsetstr)
        data['VarRef'].append(rsidstr)
    df = pandas.DataFrame(data)
    df.to_csv('marker.tsv', sep='\t', index=False)


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('--configfile', default='../../config.json')
    parser.add_argument('data', help='Table 1 from Voskoboinik et al 2018')
    return parser


def main(args):
    with open(args.configfile, 'r') as fh:
        config = json.load(fh)
        vcf = config['dbsnp38']
        rsidx = config['dbsnp38'].replace('.vcf.gz', '.rsidx')
    markernames = list(get_marker_names(args.data))
    regions = list(get_regions(args.data))
    for name, region in zip(markernames, regions):
        get_variants(name, region, config['1kgp_dir'])
        filter_variants(name)
        transfer_variants(name, vcf, rsidx)
    finalize_marker_definitions(markernames)


if __name__ == '__main__':
    main(cli().parse_args())
