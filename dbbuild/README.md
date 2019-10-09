# MicroHapDB Build

MicroHapDB is a collection of microhaplotype marker and allele frequency data integrated from a variety of public sources.
Most users will be primarily interested in MicroHapDB's data tables and access methods (CLI and Python API).
However, a bit of preliminary work was required to prepare each independent data source for compilation into a single database.
These details are described here, and those who are especially curious or who wish to reproduce this work should find everything they need here.


## Sources

Each source of microhap data has a dedicated directory located in `dbbuild/sources/`.
As of this writing, this includes `dbbuild/sources/alfred2018/`, `dbbuild/sources/vandergaag2018/`, and `dbbuild/sources/staadig2019/`.
Re-compiling the database to include additional data sources is a matter of creating a new directory<sup>1</sup> in `dbbuild/sources/`, populating it with appropriately formatted files (described in section "**Data Required for Each Source**"), and executing the database build workflow (described in section "**Running the Database Build**").

<sup>1</sup>Directories should have descriptive names, corresponding to a database name or the first author of a corresponding publication, followed by the year of access or publication.


## Data Required for Each Source

Each source directory must contain 4 data files.
A description of each required file is given below.
Raw original data files and code to transform the data into the required format should also be included.

- `marker.tsv`
- `population.tsv`
- `frequency.tsv`
- `source.txt`

### `marker.tsv`

The `marker.tsv` file contains the microhaplotype marker definition.
This includes the following fields.

- `Name`: a symbol following the `mh<chromosome><initials>-<identifier>` nomenclature proposed in [Kidd (2016)](https://dx.doi.org/10.1186/s40246-016-0078-y)
- `Reference`: the version of the human reference genome used as a variant coordinate system; at the moment, only GRCh38 is supported
- `Chrom`: the chromosome on which the marker is found
- `Offsets`: offsets (from the first nucleotide on the chromosome) of the variants that define the marker, separated by commas; note that the first nucleotide of the chromosome has an offset of `0`, the second nucleotide has an offset of `1`, and so on
- `Xref`: an identifier (or comma-separated list of identifiers) used to cross-reference the marker in other databases; this field is optional

For example, the first few lines of the `marker.tsv` for ALFRED looks like this.

```
Name	Reference	Chrom	Offsets	Xref
mh01KK-172	GRCh38	chr1	1551453,1551522,1551678	SI664721A
mh01KK-172	GRCh38	chr1	3826567,3826754,3826785,3826826	SI664246C
mh01KK-172	GRCh38	chr1	4167403,4167500,4167563,4167573	SI664548H
```

### `population.tsv`

The `population.tsv` file contains a description of all the populations for which microhap allele frequency data is available.
It includes the following fields.

- `Name`: a free-text name, intended to be human readable
- `Xref`: an identifier (or comma-separated list of identifiers) used to cross-reference the population in other databases; at least one Xref identifier is required

For example, the first few lines of the `population.tsv` for ALFRED look like this.

```
Name	Xref
Adygei	SA000017I
African Americans	SA000101C
African Americans	SA004047P
Afro-Caribbeans	SA004242M
```

**NOTE**: population names need not be unique, but `Xref` identifiers must be unique across all sources.

### `frequency.tsv`

The `frequency.tsv` file contains population allele frequency data for the microhaplotype markers.
It includes the following fields.

- `Marker`: the `Name` of the marker
- `Population`: a `Xref` identifer of the population
- `Allele`: the allele of each variant in the microhap, separated by dashes
- `Frequency`: the frequency of the allele in the specified population (a real number between 0.0 and 1.0)

For example, the first few lines of the `frequency.tsv` for ALFRED look like this.

```
Marker	Population	Allele	Frequency
mh17KK-014	SA000001B	G,C,T	0.01
mh17KK-014	SA000001B	G,C,C	0.0
mh17KK-014	SA000001B	G,A,C	0.0
mh17KK-014	SA000001B	C,C,T	0.0
mh17KK-014	SA000001B	C,C,C	0.99
mh17KK-014	SA000002C	C,A,C	0.0
```

### `source.txt`

This file should contain a single unique identifier corresponding to the source of the data.
If the data is associated with a single publication, for example, the Digital Object Identifier (DOI) of that publication is the most appropriate identifier.
If a single DOI is not a suitable identifier, an alternative label may be used.
For example:

- The allele frequency data published on the ALFRED database in December, 2018 does not correspond to a single publication.
  There is substantial overlap with the data described in [Kidd et al. (2018)](https://doi.org/10.1002/elps.201800092), but the latest ALFRED data includes additional markers and populations.
  MicroHapDB simply uses the label `ALFRED` for this source of data.
- A poster presented at ISFG 2019 described an evaluation of 45 microhaps on a Swedish population.
  The authors were kind enough to share their data and gave permission have it included in MicroHapDB.
  Since there is not yet any publication associated with the data, MicroHapDB uses the label `ISFG2019:P597` for this source of data, (poster 597 at ISFG 2019).


## Running the Database Build

More soon!


----------

# Appendix


## Prerequisites

Building MicroHapDB from scratch requires the following software and data.

- Python 3 (tested with 3.6, but will probably work with earlier 3.x versions)
- Pandas (tested with 0.23.4)
- Snakemake (tested with 5.1.5)
- dbSNP database (VCF format, gzip compressed)
- tabix (tested with 1.9)

For convenience, let's assume for the rest of this manual that the path to the dbSNP database is stored in the environmental variable `dbsnp`.

```
export dbsnp=/path/to/the/file/dbSNP_GRCh38.vcf.gz
```

The build procedure has been tested on the Linux and Mac OS X operating systems.


## Data download

> **NOTE**: This step should not need to be repeated and is described only for the purpose of full disclosure.

Important variant information for ALFRED microhaplotypes are unavailable in convenient summary form.
The `ALFRED.Snakefile` workflow was used to retrieve HTML summary pages for microhap loci over the network from the ALFRED database.

```
snakemake --snakefile ALFRED.Snakefile -p loci
```

These files are stored in the `alfred/downloads/locus-detail/` directory for scraping by the main data processing workflow.
All other data required for the build has been downloaded as described in the `alfred/` and `lovd/` directories.

> **ANOTHER NOTE**: As of December 2018, funding for the ALFRED database is set to expire.
> It is uncertain how long this build step will work.
> Fortunately the data have been captured and stored here to enable future reference.


## Data processing

The main build procedure is implemented in `Snakefile` and scrapes, cleans, formats, and cross-references data from ALFRED, LOVD, and dbSNP.
It is invoked like so.

```
snakemake tables --config dbsnp=$dbsnp
```

To speed up the build with multiple processes:

```
snakemake tables --cores 8 --config dbsnp=$dbsnp
```


## Developer note

Early versions of the software used the term `locus` liberally.
This has since been replaced by `marker` in the main code base, but `locus` and `loci` are still pervasive throughout the DB build code.
The two terms can be used interchangeably.
