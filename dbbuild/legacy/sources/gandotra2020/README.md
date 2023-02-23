# mMHseq Panel

## Citations

Gandotra N, Speed WC, Qin W, Tang Y, Pakstis AJ, Kidd KK, Scharfe C (2020) Validation of novel forensic DNA markers using multiplex microhaplotype sequencing. *Forensic Science International: Genetics*, **47**:102275, [doi:10.1016/j.fsigen.2020.102275](https://doi.org/10.1016/j.fsigen.2020.102275).


## Build process

Run the following command from the `dbbuild/sources/gandotra2020/` directory to compile the data into the table format required by MicroHapDB.

```bash
# Compile markers
./vcfadd.py manual.tsv databases/dbSNP/dbSNP_GRCh37.{vcf.gz,rsidx} databases/dbSNP/dbSNP_GRCh38.{vcf.gz,rsidx} auto.tsv novel.grch37.bed
liftover auto.tsv hg19ToHg38.over.chain.gz novel.grch38.bed novel-unmapped.bed
./novelcoords.py auto.tsv novel.grch37.bed novel.grch38.bed complete.tsv
./crosschecks.py complete.tsv databases/hg37.fa databases/hg38.fa
./prep.py complete.tsv marker.tsv

# Compile frequencies
./compile_freqs.py original/ > frequency.tsv
```


## Appendix

### Pre-processing


The `manual.tsv` file contains marker data collected manually by two individuals.
Preliminary manual results were cross-compared with each other as well as to computationally inferred results to identify and resolve discrepancies.

Frequency data is unavailable for bulk download from the mMHseq Shiny app.
The web driver script `mMHseq_download.py` was written to automatically collect all CSV frequency data from the web site.
The `mMHseq_check.py` script is available to check the integrity of the downloads.
These scripts were used to populate the CSV files in the `original/` directory
