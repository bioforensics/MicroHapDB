#!/usr/bin/env bash

for vcf in ALL.chr*.vcf.gz; do
    prefix=$(basename $vcf .vcf.gz)
    rsidx=${prefix}.rsidx
    rsidx index $vcf $rsidx
done
