# Large-scale Identification Panel

## Citations

de la Puente M, Phillips C, Xavier C, Amigo J, Carracedo A, Parson W, Lareu MV (2020) Building a custom large-scale panel of novel microhaplotypes for forensicidentification using MiSeq and Ion S5 massively parallel sequencing systems. *FSI: Genetics*, 45:102213, [doi:10.1016/j.fsigen.2019.102213](https://doi.org/10.1016/j.fsigen.2019.102213).

## Build Process

The following software is required to compile the published data into the table format required by MicroHapDB.

- Python 3
- [Snakemake][]

To build the TSV files required by MicroHapDB, run the following command from the `dbbuild/sources/vandergaag2018/` directory.

```
snakemake -p all
```

## Manual Pre-processing / known issues

???

nors --> rs772115763
rs74898010 --> rs73151289
rs28970291 --> rs4076758
rs72629020 --> rs36190610

[Snakemake]: https://snakemake.readthedocs.io/en/stable/
