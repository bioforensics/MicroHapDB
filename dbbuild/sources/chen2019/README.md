# Microhap panel for mixtures (Chen 2019)

## Citations

Chen P, Deng C, Li Z, Pu Y, Yang J, Yu Y, Li K, Li D, Liang W, Zhang L, Chen F (2019) A microhaplotypes panel for massively parallel sequencing analysis of DNA mixtures. *FSI: Genetics*, 40:140-149, [doi:10.1016/j.fsigen.2019.02.018](https://doi.org/10.1016/j.fsigen.2019.02.018).

## Build Process

The following software is required to compile the published data into the table format required by MicroHapDB.

- Python 3
- [Snakemake][]

Run the following command from the `dbbuild/sources/chen2019/` directory to compile the data into the table format required by MicroHapDB.

```bash
snakemake --config dbsnp=/path/to/dbSNP.vcf.gz rsidx=/path/to/dbSNP.rsidx -p all
```

## Manual Pre-processing

The file  `prelim-marker-variants.xslx` is a manually reformatted version of Table 1 from the manuscript.
The file `prelim-allele-freq.xlsx` is a manually reformatted version of the allele frequency table included in the supplementary data (`mmc6.xlsx`).
Credit goes to Rebecca Mitchell for this work!

These files were manually exported to TSV format to facilitate subsequent processing.


## Problematic allele

The allele frequency table included in the supplement includes an allele designated simply as "6" for the marker MH06CP003.
This marker is not one of the 11 markers described in the paper, and had been recorded previously in the ALFRED database—see https://alfred.med.yale.edu/alfred/recordinfo.asp?UNID=SI664884K.
The only allele for this marker present in ALFRED but missing from the frequency table is GAC, and we have thus substituted this for "6" in the manually reformatted allele frequency table.


[Snakemake]: https://snakemake.readthedocs.io/en/stable/
[rsidx]: https://github.com/bioforensics/rsidx
