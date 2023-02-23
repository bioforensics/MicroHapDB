# MicroHapDB Database Construction

The contents of the MicroHapDB database were integrated from a variety of public sources.
Considerable care was taken to capture original data and document the process of reformatting and streamlining this data into a format that could be integrated into MicroHapDB.
This directory describes how the core MicroHapDB tables were constructed, and provides a model for integrating marker/population/frequency data from additional sources.

> Since its inception, MicroHapDB has been envisioned as a resource requiring minimal infrastructure.
> The database interface does not require a central server, but is installed as part of a software package run locally.
> The contents of the database are stored in simple easy-to-use tables using the popular CSV format.
> The database construction process is automated and documented, enabling any reasonably competent bioinformatics practitioner to repopulate the database from scratch using primary data sources and commodity computing hardware.
> At some point in the future, MicroHapDB may support a Web-based interface, but at the moment this is a lower priority.

If MicroHapDB is missing important data, a request to include this data can be submitted using the [MicroHapDB issue tracker](https://github.com/bioforensics/MicroHapDB/issues/new).
Contributions from the community in the form of GitHub pull requests are also welcome!


## Database Sources

Microhap marker, population, and frequency data is organized by source for the database build.
Here, "source" refers to a paper describing a specific published study.
Each source has a dedicated directory located in `dbbuild/sources/.`

Adding new data to MicroHapDB requires adding a new "source" directory, populating it with data files in the expected format (described below), and executing the database build workflow.
By convention, the name of this directory is the surname of the first author of the study followed by the publication year (e.g. `kidd2018`, `delapuente2020`, `wu2021`).

Each data source must include a `source.json` file with metadata describing the study.
Additionally, each source must include at least one of the following files, but can include as many appropriate depending on available data.

- `marker.csv`
- `population.csv`
- `frequency.csv`
- `indels.csv` (required only if one ore more marker definitions include indel variants)

If these files are populated manually, no other files are required.
However, for some sources, MicroHapDB includes code for converting data from its original *ad hoc* format to the format expected by MicroHapDB.

### `source.json`

This file contains very basic metadata for the study.
It includes the following fields.

- `name`: a unique name for the data source; by convention, matches the directory name
- `year`: publication year
- `doi`: DOI of the corresponding paper
- `description`: free-text description of the data source

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

### `population.csv`

The `population.csv` file contains a description of any population groups for which MH allele frequency data is available from this study.
This includes the following fields.

- `ID`: a unique identifier for this population across all sources
- `Name`: a free-text description of the population, intended to be human readable
- `Xref`: optional cross-reference

For example, the first few lines of the `population.tsv` for the original 1000 Genomes Project data (Auton 2015) look like this.

```csv
ID,Name,Xref
CHB,"Han Chinese in Beijing, China",SA004058R
JPT,"Japanese in Tokyo, Japan",SA004060K
CHS,Southern Han Chinese,SA004059S
CDX,"Chinese Dai in Xishuangbanna, China",SA004238R
KHV,"Kinh in Ho Chi Minh City, Vietnam",SA004249T
CEU,Utah Residents (CEPH) with Northern and Western European Ancestry,SA004250L
TSI,Toscani in Italia,SA004057Q
FIN,Finnish in Finland,SA004049R
GBR,British in England and Scotland,SA004050J
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

### `indel.csv`

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

> **NOTE**: MicroHapDB uses this indel metadata to display MH sequences correctly.
> Population frequency estimates for microhaps with indels may not be accurate for all data sources, and it is not recommended to rely on frequency estimates for these markers.


## Rebuilding the Database from Scratch

With software dependencies and databases in place, the database can be rebuilt with the following command (run within the `dbbuild/` directory).

```
./build.py | tee build-summary.txt
```

The CSV files created during the build should then be copied to the `microhapdb/data/` directory.

### Dependencies

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

### Databases

- dbSNP
    - .vcf.gz, .vcf.gz.tbi, and .rsidx files
    - GRCh37 and GRCh38
- 1000 Genomes VCFs
- UCSC liftover chain files
    - hg19ToHg38
    - hg38ToHg19
