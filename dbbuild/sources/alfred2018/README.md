# Allele Frequency Database (ALFRED)

## Citations

Rajeevan H, Soundararajan U, Kidd JR, Pakstis AJ, Kidd KK (2012) ALFRED: an allele frequency resource for research and teaching. *Nucleic Acids Research*, 40(D1): D1010-D1015, [doi:10.1093/nar/gkr924](https://doi.org/10.1093/nar/gkr924).

Kidd KK, Pakstis AJ, Speed WC, Lagace R, Wootton S, Chang J (2018) Selecting microhaplotypes optimized for different purposes. *Electrophoresis*, [doi:10.1002/elps.201800092](https://doi.org/10.1002/elps.201800092).

Kidd KK, Rajeevan H (2018) ALFRED data download. *The Allele Frequency Database*, https://alfred.med.yale.edu/alfred/selectDownload/Microhap_alleleF_198.txt. Accessed December 7, 2018.

## Build Process

### Step 1: download ALFRED marker and frequency data (OPTIONAL)

Microhap allele frequencies are available from ALFRED in a single text file.
However, the marker definitions are spread across many HTML pages.
These HTML files have been downloaded to `dbbuild/sources/alfred2018/marker-detail/`.
**Downloading them again shouldn't be necessary**, but for the sake of the curious the download process can be repeated by running the following command from the `dbbuild/sources/alfred2018/` directory.

```
snakemake -s download.Snakefile frequencies markers
```

### Step 2: compile data tables

Run the following command from the `dbbuild/sources/alfred2018/` directory to compile the data into the table format required by MicroHapDB.

```
snakemake --configfile ../../config.json -p all
```


## Appendix

### Known Issues

The ALFRED database publishes 1GKP-derived frequencies for some markers.
For many markers there is perfect agreement between ALFRED and MicroHapDB, but in some cases there are slight or notable differences.
After correspondence with the ALFRED curators, we suspect that these differences are due to their use of PHASE over all aggregated data sets to statistically infer haplotypes, while MicroHapDB relies entirely on haplotypes as published in the 1KGP data.

1KGP-derived frequencies published by ALFRED are stored here in `frequency-1kgp.tsv` but are not integrated into the final MicroHapDB database.
