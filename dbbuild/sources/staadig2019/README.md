# Linköping Panel Evaluation (Sweden)

## Citations

Staadig A, Tillmar A (2019) Evaluation of microhaplotypes—A promising new type of forensic marker. *The 28th Congress of the International Society for Forensic Genetics*, P597.

## Build Process

Run the following command from the `dbbuild/sources/alfred2018/` directory to compile the data into the table format required by MicroHapDB.

```
snakemake --configfile ../../config.json -p all
```


## Appendix

### Manual Pre-processing

The files `marker-variants.tsv` and `allele-frequencies.tsv` were exported to TSV from the corresponding `Position` and `Frequencies` tabs of the original data file furnished by the authors (in `original/Microhaplotype positions and frequencies.xlsx`).

The population ID was created by appending the output of `echo $'ISFG2019:P597\tSwedish' | md5 | cut -c 1-10` to the prefix `MHDBP-`.


[Snakemake]: https://snakemake.readthedocs.io/en/stable/
[rsidx]: https://github.com/bioforensics/rsidx
