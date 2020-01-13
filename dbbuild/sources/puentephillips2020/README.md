# Large-scale Identification Panel

## Citations

de la Puente M, Phillips C, Xavier C, Amigo J, Carracedo A, Parson W, Lareu MV (2020) Building a custom large-scale panel of novel microhaplotypes for forensicidentification using MiSeq and Ion S5 massively parallel sequencing systems. *FSI: Genetics*, 45:102213, [doi:10.1016/j.fsigen.2019.102213](https://doi.org/10.1016/j.fsigen.2019.102213).

## Build Process

The following software is required to compile the published data into the table format required by MicroHapDB.

- Python 3
- [Snakemake][]
- [rsidx][]

To build the TSV files required by MicroHapDB, run the following command from the `dbbuild/sources/puentephillips2020/` directory.

```
snakemake --config dbsnp=/path/to/dbSNP.vcf.gz rsidx=/path/to/dbSNP.rsidx -p all
```

## Manual Pre-processing

The [pdfminer](https://github.com/euske/pdfminer) package was used to extract text from the Supplementary File S1 (PDF) file.
One variant in the file was marked as "nors", but manual examination confirmed that rs772115763 refers to the SNP at the correct position with the expected alleles.

Three additional SNPs are labeled in File S1 with rsIDs that have been deprecated and merged into newer rsIDs.
The `text2table.py` script makes the following rsID replacements.

```
rs74898010 --> rs73151289
rs28970291 --> rs4076758
rs72629020 --> rs36190610
```

[Snakemake]: https://snakemake.readthedocs.io/en/stable/
[rsidx]: https://github.com/bioforensics/rsidx
