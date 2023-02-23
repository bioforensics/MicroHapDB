# USC Large-scale Identification Panel

## Citations

de la Puente M, Phillips C, Xavier C, Amigo J, Carracedo A, Parson W, Lareu MV (2020) Building a custom large-scale panel of novel microhaplotypes for forensic identification using MiSeq and Ion S5 massively parallel sequencing systems. *FSI: Genetics*, 45:102213, [doi:10.1016/j.fsigen.2019.102213](https://doi.org/10.1016/j.fsigen.2019.102213).

## Build Process

Run the following command from the `dbbuild/sources/delapuente2020/` directory to compile the data into the table format required by MicroHapDB.

```
snakemake --configfile ../../config.json --cores 1 -p all
```

## Appendix

## Manual Pre-processing

The [pdfminer](https://github.com/euske/pdfminer) package was used to extract text from the Supplementary File S1 (PDF) file.
One variant in the file was marked as "nors."
Manual examination confirms that rs772115763 refers to the SNP at the correct position with the expected alleles.
The "nors" string is substituted with rs772115763 by the `text2table.py` script.

[Snakemake]: https://snakemake.readthedocs.io/en/stable/
[rsidx]: https://github.com/bioforensics/rsidx
