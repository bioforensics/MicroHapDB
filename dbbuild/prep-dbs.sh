#!/usr/bin/env bash
set -eo pipefail

download_1kgp()
{
    local rootdir=$1
    local rooturl="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    mkdir -p ${rootdir}/1000Genomes/
    for i in {1..22}; do
        url="${rooturl}/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
        filename="${rootdir}/1000Genomes/chr${i}.vcf.gz"
        wget -O ${filename} ${url}
        wget -O ${filename}.tbi ${url}.tbi
    done
    url="${root}/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
    filename="${rootdir}/1000Genomes/chrX.vcf.gz"
    wget -O ${filename} ${url}
    wget -O ${filename}.tbi ${url}.tbi
}

index_1kgp()
{
    local rootdir=$1
    local numjobs=$2
    parallel --gnu --jobs $2 rsidx index ${rootdir}/1000Genomes/chr{}.vcf.gz ${rootdir}/1000Genomes/chr{}.rsidx ::: {1..22} X
}

download_grch38()
{
    local rootdir=$1
    wget -O ${rootdir}/hg38.fa.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    gunzip ${rootdir}/hg38.fa.gz
}

download_dbsnp()
{
    local rootdir=$1
    local rooturl="https://ftp.ncbi.nlm.nih.gov/snp/organisms/"
    mkdir -p ${rootdir}/dbSNP/
    url="${rooturl}/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz"
    filename="${rootdir}/dbSNP/dbSNP_GRCh37.vcf.gz"
    wget -O ${filename} ${url}
    wget -O ${filename}.tbi ${url}.tbi
    url="${rooturl}/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz"
    filename="${rootdir}/dbSNP/dbSNP_GRCh38.vcf.gz"
    wget -O ${filename} ${url}
    wget -O ${filename}.tbi ${url}.tbi
}

index_dbsnp()
{
    local rootdir=$1
    local numjobs=$2
    parallel --gnu --jobs $numjobs rsidx index ${rootdir}/dbSNP/dbSNP_{}.vcf.gz ${rootdir}/dbSNP/dbSNP_{}.rsidx ::: GRCh37 GRCh38
}

dbdir=${1:-databases}
cores=${2:-1}

download_grch38 $dbdir
download_1kgp $dbdir
download_dbsnp $dbdir

if [ $cores -le 4 ]; then
    index_dbsnp $dbdir $cores
    index_1kgp $dbdir $cores
else
    index_dbsnp $dbdir 2 &
    index_1kgp $dbdir $((cores - 2))
fi
