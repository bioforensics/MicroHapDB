# MicroHapDB Database Construction

## Overview

This directory describes the MicroHapDB database build procedure and provides a model for integrating marker, population, and frequency data from additional sources.
Considerable care was taken to capture original data and document the process of reformatting and streamlining this data for aggregation in MicroHapDB.
This directory describes how the core MicroHapDB tables were constructed, and provides a model for integrating marker and frequency data from additional sources.

If MicroHapDB is missing important data, a request to include this data can be submitted using the [MicroHapDB issue tracker](https://github.com/bioforensics/MicroHapDB/issues/new).
Contributions from the community in the form of GitHub pull requests are also welcome!


## Primary Data Sources

The `sources/` directory contains a dedicated subdirectory for each primary data source integrated into MicroHapDB.
Each subdirectory in turn contains "raw" primary data files as published in the literature, often in the form of tables, text files, or spreadsheets.
It also contains the code used to reformat the primary data into a consistent format usable by MicroHapDB.

Adding new data sources to MicroHapDB typically does not require re-processing any of the data sources already present.
The code and documentation for each primary data source is provided in the spirit of full disclosure, but also establishes the complete provenance of each data point, facilitates reproducibility, and allows any errors in data processing to be corrected rapidly.
The goal is that—if ever needed, heaven forbid—any reasonbly capable bioinformatician with an interest in microhaplotypes could pick up and take over stewardship of this community resource.

> Some sources include dedicated documentation in a README file, providing build instructions tailored specifically for that data source.
> For sources lacking a README, the default build process is to run `snakemake -c1` in the directory.


## Rebuilding the Database from Scratch

The MicroHapDB database can be rebuilt with the following command in the `dbbuild/` directory.

```
./build.py databases/dbSNP/ databases/chains/ | tee build-summary.txt
```

The arguments provided to the build script will depend on the location of the dbSNP files and liftover chain files on the system.
Running `./build.py --help` should provide helpful guidance.

If the build is successful, the updated data tables can be copied to the main data directory with the following command.

```
cp *.csv* ../microhapdb/data/
```

Note that if the build script is run to integrate new marker definitions, 1000 Genomes Project population frequency estimates for those new markers will not be computed without re-running the `sources/byrskabishop2022/` build procedure after updating its `marker-latest.csv` file.
The `./build.py` script must be run before the `sources/byrskabishop2022/` build to provide an up-to-date marker file, and then it must be run again after the `sources/byrskabishop2022/` build to aggregate the newly computed frequency and $A_e$ values.

See the appendices at the end of this document for details about which software dependencies and databases are required for the database build.


## Data Types

MicroHapDB's core data is comprised of marker definitions and allele frequency estimates.
Each primary data source must contain at least one of these data types, as well as source metadata.
The required format for this data is described below.

### `source.json`

This file contains very basic metadata for the study.
It includes the following fields.

- `name`: a unique name for the data source; by convention, matches the directory name
- `year`: publication year
- `doi`: DOI of the corresponding paper
- `description`: free-text description of the data source
- `order`: optionally, specify the order of precedence this source should have in its publication year; when not provided, sources can have an arbitrary precedence within a given year

For example, the `source.json` file for Kidd (2018) looks like this.

```json
{
    "name": "Kidd2018",
    "year": 2018,
    "doi": "10.1002/elps.201800092",
    "description": "198 MH loci collected from the ALFRED database, along with allele frequencies for 83 population groups"
}
```

### `marker.csv`

The `marker.csv` file contains microhaplotype marker definitions.
This includes the following fields.

- `Name`: a symbol following the `mh<chromosome><initials>-<identifier>` nomenclature proposed in [Kidd (2016)](https://dx.doi.org/10.1186/s40246-016-0078-y)
- `Xref`: an identifier (or semicolon-separated list of identifiers) used to cross-reference the marker in other databases; this field can be empty, but if it includes one or more identifiers, each of these identifiers should be unique to this marker across all sources
- `NumVars`: the number of SNPs that define this marker
- `Refr`: the human genome reference assembly version to which the following coordinates correspond; must be "GRCh37" (for hg19) or "GRCh38"; *optional if a full RSID list is provided*
- `Chrom`: the chromosome on which the marker is found
- `Positions`: a semicolon-separated list of SNP coordinates; *optional if a full RSID list is provided*
- `VarRef`: semicolon-separated list of variant identifiers (rsIDs) used to define this marker

For example, the first few lines of the `marker.csv` for de la Puente (2020) look like this.

```csv
Name,Xref,NumVars,Refr,Chrom,Positions,VarRef
mh01USC-1pA,,4,,chr1,,rs28503881;rs4648788;rs72634811;rs28689700
mh01USC-1pB,,3,,chr1,,rs3766785;rs1679910;rs1679911
mh01USC-1pC,,3,,chr1,,rs4472717;rs1776320;rs2211434
mh01USC-1pD,,3,,chr1,,rs6702428;rs12031966;rs6687440
mh01USC-1qA,,2,,chr1,,rs870880;rs870879
mh01USC-1qB,,4,,chr1,,rs4512600;rs10753020;rs4344286;rs34715816
mh01USC-1qC,,3,,chr1,,rs6702229;rs117262471;rs6665072
mh01USC-1qD,,3,,chr1,,rs12141154;rs74148721;rs56274766
mh02USC-2pA,,4,,chr2,,rs76708321;rs1686420;rs1686419;rs1686418
```

### `frequency.csv`

The `frequency.csv` file contains any population frequency data published by the study.
It includes the following fields.

- `Marker`: the `Name` of the marker
- `Population`: the unique identifer of the population
- `Allele`: the allele of each variant in the microhap, separated by pipe symbols
- `Frequency`: the frequency of the allele in the specified population (a real number between 0.0 and 1.0)

For example, the first few lines of the `frequency.tsv` for van der Gaag (2018) look like this.

```csv
Marker,Population,Allele,Frequency
mh06PK-24844,MHDBP-383d86606a,T|C|G|C|C|C|A|A|G|A,0.000
mh06PK-24844,MHDBP-936bc36f79,T|C|G|C|C|C|A|A|G|A,0.000
mh06PK-24844,MHDBP-3dab7bdd14,T|C|G|C|C|C|A|A|G|A,0.123
mh06PK-24844,MHDBP-383d86606a,T|C|G|C|C|T|A|A|G|G,0.566
mh06PK-24844,MHDBP-936bc36f79,T|C|G|C|C|T|A|A|G|G,0.586
mh06PK-24844,MHDBP-3dab7bdd14,T|C|G|C|C|T|A|A|G|G,0.425
mh06PK-24844,MHDBP-383d86606a,C|C|G|C|C|C|A|A|G|A,0.071
mh06PK-24844,MHDBP-936bc36f79,C|C|G|C|C|C|A|A|G|A,0.000
mh06PK-24844,MHDBP-3dab7bdd14,C|C|G|C|C|C|A|A|G|A,0.329
```

### `population.csv`

The `population.csv` file provides a description of any population groups for which MH allele frequency data is available from this study.
This includes the following fields.

- `ID`: a unique identifier for this population across all sources
- `Name`: a free-text description of the population, intended to be human readable
- `Xref`: optional cross-reference

For example, the first few lines of the `population.tsv` for the NYGC 1000 Genomes Project data (Byrska-Bishop 2022) look like this.

```csv
ID,Name,Xref
1KGP,1000 Genomes Aggregate,
AFR,Africa,
EUR,Europe,
SAS,South Asia,
EAS,East Asia,
CHB,"Han Chinese in Beijing, China",SA004058R
JPT,"Japanese in Tokyo, Japan",SA004060K
CHS,Southern Han Chinese,SA004059S
```

### `indel.csv`

> **NOTE**: MicroHapDB uses this indel metadata to display MH sequences correctly.
> Population frequency estimates for microhaps with indels may not be accurate for all data sources, and it is not recommended to rely on frequency estimates for these markers.

The `indel.csv` file contains insertion/deletion variant data for any markers that include indels.
It includes the following fields.

- `Marker`: the `Name` of the marker
- `VariantIndex`: the index of the variant within the microhap marker; the first variant in the marker is `0`, the next is `1`, and so on
- `Refr`: the allele represented in the reference genome
- `Alt`: the alternate allele (or comma-separated list of alternate alleles) present in the marker definition

For example, the entire `indel.csv` file for Kidd (2018) looks like this.

```
Marker,VariantIndex,Refr,Alt
mh07KK-081,0,TA,T
mh11KK-091,0,TG,T
mh22KK-064,3,AATAATT,A
```


## Appendix A: Software Dependencies

The software dependencies for the build include the following.
It should be possible to possible to perform the main build with a subset of these dependencies, but rebuilding some specific sources is only possible with all of these dependencies.
They can be installed using pip and/or conda.

- pandas
- snakemake
- pyfaidx
- rsidx
- scikit-allel
- ucsc-liftover
- selenium
- geckodriver
- intervaltree
- tqdm


## Appendix B: Required Auxiliary Data files

- dbSNP
    - .vcf.gz, .vcf.gz.tbi, and .rsidx files for both GRCh37 and GRCh38
    - info on merged records
- UCSC liftover chain files
    - hg19ToHg38
    - hg38ToHg19

The following command will download data files required for the database build.

```
snakemake -c1 -p -s download.smk -d databases/
```

Following a successful run, the command `./build.py --check` can be used to verify that the files were downloaded correctly to the expected location.

A dbSNP rsidx index must also be built for both GRCh37 and GRCh38 coordinates. Note that the following commands require many hours of run time.

```
rsidx index databases/dbSNP/dbSNP_GRCh37.vcf.gz databases/dbSNP/dbSNP_GRCh37.rsidx
rsidx index databases/dbSNP/dbSNP_GRCh38.vcf.gz databases/dbSNP/dbSNP_GRCh38.rsidx
```
