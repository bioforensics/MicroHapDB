# Linköping Panel Evaluation (Sweden)

## Citations

Staadig A, Tillmar A (2019) Evaluation of microhaplotypes—A promising new type of forensic marker. *The 28th Congress of the International Society for Forensic Genetics*, P597.

## Build Process

The following software is required to compile the published data into the table format required by MicroHapDB.

- Python 3
- [Snakemake][]
- [rsidx][]

A recent version of the dbSNP database is also required, both the VCF file and the corresponding tabix index.

Run the following commands from the `dbbuild/sources/staadig2019/` directory to compile the data into the table format required by MicroHapDB.

```bash
# Build rsidx index if it doesn't yet exist; requires > 1 hour
rsidx index /path/to/dbSNP.vcf.gz /path/to/dbSNP.rsidx

# Compile the data tables
snakemake --config dbsnp=/path/to/dbSNP.vcf.gz rsidx=/path/to/dbSNP.rsidx -p all
```

## Manual Pre-processing

The files `marker-variants.tsv` and `allele-frequencies.tsv` were exported to TSV from the corresponding `Position` and `Frequencies` tabs of the original data file furnished by the authors (in `original/Microhaplotype positions and frequencies.xlsx`).
