# 1000 Genomes Project

## Citations

Auton A, Abecasis G, Altshuler D, et al. (2015) A global reference for human genetic variation. *Nature* **526**:68–74,  [doi:10.1038/nature15393](https://doi.org/10.1038/nature15393).


## Build Process

The following software is required to compile the published data into the table format required by MicroHapDB.

- Python 3
- wget
- tabix (part of the [samtools/HTSlib project](https://github.com/samtools/htslib))
- [rsidx][]

Run the following commands from the `dbbuild/sources/1kgp/` directory to compile the data into the table format required by MicroHapDB.

```bash
# Download 1000 Genomes Project Phase 3 data
# Note: this script is also used to build the voskoboinik2018 source; theis data only needs to be downloaded once
./download.sh

# Build rsidx index for 1KGP VCFs
for vcf in ALL.chr*.vcf.gz; do
    prefix=$(basename $vcf .vcf.gz)
    rsidx=${prefix}.rsidx
    rsidx index $vcf $rsidx
done

# Get mapping of marker IDs to rsIDs
./get-marker-rsids.py > marker-rsids-MicroHapDB-0.5.tsv

# Compute allele frequency estimates
./mhfreqs.py sample-pops.tsv marker-rsids-MicroHapDB-0.5.tsv > frequency.tsv
```


## Manual Pre-processing

The `population.tsv` file was copied, pasted, and edited manually from the ISGR at https://www.internationalgenome.org/faq/which-populations-are-part-your-study/.

The `sample-pops.tsv` file was created by downloading `original/20130606_sample_info.xlsx` from [this page](https://www.internationalgenome.org/faq/which-samples-are-you-sequencing/) (URL: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx) and manually editing the relevant columns.


## Known Issues

The ALFRED database publishes 1GKP-derived frequencies for some markers.
For many markers there is perfect agreement between ALFRED and MicroHapDB, but in some cases there are slight or notable differences.
After correspondence with the ALFRED curators, we suspect that these differences are due to their use of PHASE over all aggregated data sets to statistically infer haplotypes, while MicroHapDB relies entirely on haplotypes as published in the 1KGP data.


[rsidx]: https://github.com/bioforensics/rsidx
