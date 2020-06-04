# 1000 Genomes Project

## Citations

Auton A, Abecasis G, Altshuler D, et al. (2015) A global reference for human genetic variation. *Nature* **526**:68–74, [doi:10.1038/nature15393](https://doi.org/10.1038/nature15393).


## Build Process

The following software is required to compile the published data into the table format required by MicroHapDB.

- Python 3
- wget
- tabix (part of the [samtools/HTSlib project](https://github.com/samtools/htslib))
- infocalc (included in this directory)
- [pandas][]
- [rsidx][]
- [scikit-allel][]
- [Snakemake][]

Run the following commands from the `dbbuild/sources/1kgp/` directory to compile the data into the table format required by MicroHapDB.

```bash
snakemake -p all
```

## Manual Pre-processing

The `population.tsv` file was copied, pasted, and edited manually from the ISGR at https://www.internationalgenome.org/faq/which-populations-are-part-your-study/.

The `sample-pops.tsv` file was created by downloading `original/20130606_sample_info.xlsx` from [this page](https://www.internationalgenome.org/faq/which-samples-are-you-sequencing/) (URL: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx) and manually editing the relevant columns.


## Known Issues

The ALFRED database publishes 1GKP-derived frequencies for some markers.
For many markers there is perfect agreement between ALFRED and MicroHapDB, but in some cases there are slight or notable differences.
After correspondence with the ALFRED curators, we suspect that these differences are due to their use of PHASE over all aggregated data sets to statistically infer haplotypes, while MicroHapDB relies entirely on haplotypes as published in the 1KGP data.

Some markers are defined by rsIDs that are present in recent versions of dbSNP, but have no corresponding entries in the 1KGP VCFs.
These markers are excluded from MicroHapDB's 1KGP frequency estimates. 

Some markers are defined by rsIDs that have been updated since the 1KGP Phase 3 VCFs were published.
For the purposes of estimating microhap frequencies, these variants are reverted to their old rsID.
See `get-marker-rsids.py` for details.


[pandas]: https://pandas.pydata.org
[rsidx]: https://github.com/bioforensics/rsidx
[scikit-allel]: https://scikit-allel.readthedocs.io
[Snakemake]: https://snakemake.readthedocs.io/en/stable/
