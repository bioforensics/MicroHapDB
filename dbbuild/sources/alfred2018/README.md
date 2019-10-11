# Allele Frequency Database (ALFRED)

## Citations

Rajeevan H, Soundararajan U, Kidd JR, Pakstis AJ, Kidd KK (2012) ALFRED: an allele frequency resource for research and teaching. *Nucleic Acids Research*, 40(D1): D1010-D1015, [doi:10.1093/nar/gkr924](https://doi.org/10.1093/nar/gkr924).

Kidd KK, Pakstis AJ, Speed WC, Lagace R, Wootton S, Chang J (2018) Selecting microhaplotypes optimized for different purposes. *Electrophoresis*, [doi:10.1002/elps.201800092](https://doi.org/10.1002/elps.201800092).

Kidd KK, Rajeevan H (2018) ALFRED data download. *The Allele Frequency Database*, https://alfred.med.yale.edu/alfred/selectDownload/Microhap_alleleF_198.txt. Accessed December 7, 2018.

## Build Process

The following software is required to compile the published data into the table format required by MicroHapDB.

- Python 3
- [Snakemake][]
- [rsidx][]

A recent version of the dbSNP database is also required, both the VCF file and the corresponding tabix index.

### Step 0: index dbSNP for rsid searches

Run the following command to index the dbSNP VCF file for rsid searches.
This usually requires more than an hour.

```
rsidx index /path/to/dbSNP.vcf.gz /path/to/dbSNP.rsidx
```

### Step 1: download ALFRED marker and frequency data

Microhap allele frequencies are available from ALFRED in a single text file.
However, the marker definitions are spread across many HTML pages.
These HTML files have been downloaded to `dbbuild/sources/alfred2018/marker-detail/`.
**Downloading them again shouldn't be necessary**, but for the sake of the curious the download process can be repeated by running the following command from the `dbbuild/sources/alfred2018/` directory.

```
snakemake -s download.Snakefile frequencies
snakemake -s download.Snakefile markers
```


To build the TSV files required by MicroHapDB, run the following command from the `dbbuild/sources/vandergaag2018/` directory.

```
snakemake --config dbsnp=/path/to/your/dbSNP.vcf.gz -p all
```

## Manual Pre-processing

The file `figure-S1.txt` was manually created using the supplementary PDF file from the paper (`original/mmc1.pdf`).
Genomic coordinates of each marker were determined by copying the sequence from the PDF, editing in some cases, and pasting into a UCSC BLAT search.
The offset of each variant site from the first nucleotide in the marker was manually checked, and random spot checking was used to verify the accuracy.

The file `table-S2.txt` was created by copying the table from the supplementary XLSX file from the paper (`original/mmc4.xlsx`) and pasting into a text editor.

The file `population.tsv` was created manually.


[Snakemake]: https://snakemake.readthedocs.io/en/stable/
[rsidx]: https://github.com/bioforensics/rsidx
