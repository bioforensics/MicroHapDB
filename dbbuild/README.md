# MicroHapDB Build

MicroHapDB is a collection of microhaplotype marker and allele frequency data integrated from a variety of public sources.
Most users will be primarily interested in MicroHapDB's data tables and access methods (CLI and Python API).
However, a bit of preliminary work was required to prepare each independent data source for compilation into a single database.
These details are described here, and those who are especially curious or who wish to reproduce this work should find everything they need here.


## Rebuilding the Database From Scratch: The Short Version

You'll need the following software:

- Python 3
- [Pandas][]
- [Snakemake][]
- [pyfaidx][]

And the following data:

- Human reference genome (version GRCh38; http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)

With these dependencies installed, invoke the following command in the `dbbuild/` directory to re-build the database from scratch.

```
snakemake --config refr=/path/to/hg38.fasta -p tables
```

FIXME

- [ ] finish documentation
- [ ] update testsc

FIXME


## The Long Version: Sources

Each source of microhap data has a dedicated directory located in `dbbuild/sources/`.
As of this writing, this includes `dbbuild/sources/alfred2018/`, `dbbuild/sources/vandergaag2018/`, and `dbbuild/sources/staadig2019/`.
Re-compiling the database to include additional data sources is a matter of creating a new directory<sup>1</sup> in `dbbuild/sources/`, populating it with appropriately formatted files (described in section "**Data Required for Each Source**"), and executing the database build workflow (described in section "**Running the Database Build**").

<sup>1</sup>Directories should have descriptive names, corresponding to a database name or the first author of a corresponding publication, followed by the year of access or publication.


## Data Required for Each Source

Each source directory must contain 4-5 data files.
A description of each required file is given below.

- `marker.tsv`
- `population.tsv`
- `frequency.tsv`
- `source.txt`
- `indels.tsv` (required only if marker definitions include insertion/deletion variants)

In addition to these 4-5 files, it is expected that raw original data files will be provided, along with code to transform the data into the required format and instructions for running the code.
No constraints are imposed regarding how these supporting files are named, what programming language is used, etc.—only that the pre-processing procedure can be replicated without requiring an unreasonable amount of time or effort.

### `marker.tsv`

The `marker.tsv` file contains the microhaplotype marker definition.
This includes the following fields.

- `Name`: a symbol following the `mh<chromosome><initials>-<identifier>` nomenclature proposed in [Kidd (2016)](https://dx.doi.org/10.1186/s40246-016-0078-y)
- `Xref`: an identifier (or comma-separated list of identifiers) used to cross-reference the marker in other databases; this field can be empty, but if it includes one or more identifiers, each of these identifiers should be unique to this marker across all sources
- `Reference`: the version of the human reference genome used as a variant coordinate system; at the moment, only `GRCh38` is supported
- `Chrom`: the chromosome on which the marker is found
- `Offsets`: offsets (from the first nucleotide on the chromosome) of the variants that define the marker, separated by commas; note that the first nucleotide of the chromosome has an offset of `0`, the second nucleotide has an offset of `1`, and so on
- `VarRef`: a list of variant identifiers (such as rsIDs) associated with this marker; optional

For example, the first few lines of the `marker.tsv` for ALFRED looks like this.

```
Name	Xref	Reference	Chrom	Offsets	VarRef
mh17KK-014	SI664726F	GRCh38	chr17	4497060,4497088,4497096	rs333113,rs8074965,rs11657785
mh03KK-150	SI664563E	GRCh38	chr3	131927127,131927156,131927242,131927311	rs1225051,rs1225050,rs1225049,rs1225048
mh04KK-010	SI664564F	GRCh38	chr4	1985210,1985244	rs3135123,rs495367
mh04KK-011	SI664565G	GRCh38	chr4	37857268,37857332	rs6855439,rs6531591
```


### `population.tsv`

The `population.tsv` file contains a description of all the populations for which microhap allele frequency data is available.
It includes the following fields.

- `ID`: a unique identifier for this population across all sources
- `Name`: a free-text description of the population, intended to be human readable

For example, the first few lines of the `population.tsv` for ALFRED look like this.

```
ID	Name
SA000001B	Han
SA000002C	Ami
SA000003D	Hakka
SA000005F	Biaka
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
SI664726F	SA000001B	G,C,T	0.01
SI664726F	SA000001B	G,C,C	0.0
SI664726F	SA000001B	G,A,C	0.0
SI664726F	SA000001B	C,C,T	0.0
SI664726F	SA000001B	C,C,C	0.99
SI664726F	SA000001B	C,A,C	0.0
SI664726F	SA000002C	G,C,T	0.0
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


### `indel.tsv`

The `indel.tsv` file contains insertion/deletion variant data for any markers that include indels.
It includes the following fields.

- `Marker`: the `Name` of the marker
- `VariantIndex`: the index of the variant within the microhap marker; the first variant in the marker is `0`, the next is `1`, and so on
- `Refr`: the allele represented in the reference genome
- `Alt`: the alternate allele (or comma-separated list of alternate alleles) present in the marker definition

For example, the `indel.tsv` file for ALFRED looks like this.

```
Marker	VariantIndex	Refr	Alt
mh07KK-081	0	TA	T
mh11KK-091	0	TG	T
mh22KK-064	3	AATAATT	A
```


[Pandas]: https://pandas.pydata.org
[Snakemake]: https://snakemake.readthedocs.io/en/stable/
[pyfaidx]: https://github.com/mdshw5/pyfaidx
