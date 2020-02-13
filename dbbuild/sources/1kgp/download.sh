#!/usr/bin/env bash
set -e

root="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"

for i in {1..22}; do
    filename="${root}/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    wget ${filename}.tbi
    wget ${filename}
done

filename="${root}/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
wget ${filename}.tbi
wget ${filename}
