# Downloads for Kidd 2018 / ALFRED

The file `Microhap_alleleF_198.txt` contains population frequency estimates for 198 MH markers in 96 global populations.
The files in the `marker-detail/` directory contain the full marker definition as published in the ALFRED database.
These files were most recently downloaded in November 2022, but have been static since 2019 or earlier.

All required files should be present, but the following commands—run in this directory—could be used to download them again from scratch.

```
snakemake
snakemake --cores 4
```

The first `snakemake` command will download the allele frequency file.
The second command will read marker info from this file to download all of the marker definition files.
