# Short hypervariable microhaplotypes (van der Gaag **et al.** 2018)

## Citations

van der Gaag KJ, de Leeuw RH, Laros JFJ, den Dunnen JT, de Knijff P (2018) Short hypervariable microhaplotypes: A novel set of very short high discriminating power loci without stutter artefacts. *Forensic Science International: Genetics*, 35:169-175, [doi:10.1016/j.fsigen.2018.05.008](https://doi.org/10.1016/j.fsigen.2018.05.008).

## Build Process

The following software is required to compile the published data into the table format required by MicroHapDB.

- Python 3
- [Snakemake][]
- tabix (from [htslib][])

A recent version of the dbSNP database is also required, both the VCF file and the corresponding tabix index.

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
[htslib]: https://github.com/samtools/htslib