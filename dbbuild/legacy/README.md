# Extending MicroHapDB or Rebuilding from Scratch

MicroHapDB is a collection of microhaplotype marker and allele frequency data integrated from a variety of public sources.
Folks interested only in using public data will likely be most interested in the data files in the `microhapdb/data/` directory, and the Python wrapper MicroHapDB provides for querying these files.
This directory describes how those data files were constructed, and includes instructions for rebuilding the database from scratch
We expect this will be helpful for folks who want to reproduce the database build or who want to integrate additional public or private data.

If you would like to rebuild your personal copy of MicroHapDB to **include private microhap data**, see the documentation below.
Rebuilding your personal copy of MicroHapDB does not alter the public MicroHapDB database.

If you would like to **submit new data to be included in MicroHapDB publicly**, let us know by opening a new thread on [MicroHapDB's issue tracker](https://github.com/bioforensics/MicroHapDB/issues/new) and attaching your data files.
Alternatively, you can submit a pull request to the MicroHapDB Github repository.


## Rebuilding the Database From Scratch: The Short Version

Building MicroHapDB from scratch requires installing several software packages.
Conda provides the most convenient way to install these.

```
conda create -c bioconda --name microhapdb -y \
    python=3.7 pandas snakemake pyfaidx rsidx scikit-allel \
    ucsc-liftover selenium geckodriver
conda activate microhapdb
```

Next, building MicroHapDB from scratch also requires the human reference genome, dbSNP, and the 1000 Genomes Project Phase 3 data.
These can be downloaded and indexed using with the `prep-dbs.sh` script.
**Note**: this can take several hours, depending on the speed of the Internet connection, computer processors, and other factors.

```
git clone https://github.com/bioforensics/MicroHapDB.git  # If you haven't already done so
cd MicroHapDB/dbbuild/
./prep-dbs.sh
```

Finally, with these data sets in place, MicroHapDB can be built using `Snakemake`.

```
snakemake --configfile config.json --cores 1 -p tables
```

If the Snakemake build process completes successfully, copy the newly created database tables to MicroHapDB's main data directory to complete the update!

```
cp *.tsv ../microhapdb/data/
```

> **NOTE**: If you run the Snakemake command and it says there is nothing to be done, you may need to delete some of the .tsv files to trigger a rebuild.

If adding new markers to MicroHapDB, you'll want to perform a second database rebuild.
The first build and copy adds the markers to the database.
You'll then want to go to the `dbbuild/sources/1kgp/` directory, recompute A<sub>e</sub> and I<sub>n</sub> statistics, and then run the root build process again.

```
cd sources/1kgp/
rm frequency.tsv marker-informativeness.tsv marker-aes.tsv marker-fst.tsv marker-rsids-MicroHapDB-latest.tsv
snakemake --cores 1 -p all
cd ../../
rm marker.tsv frequency.tsv
snakemake --configfile config.json --cores 1 -p tables
cp *.tsv sources/1kgp/marker-aes.tsv ../microhapdb/data/
```

> **NOTE**: The `marker.tsv` file contains A<sub>e</sub> values aggregated over 26 global populations, whereas the `sources/1kgp/marker-aes.tsv` contains per-population A<sub>e</sub> values for each marker.


## The Long Version: Sources

Each source of microhap data has a dedicated directory located in `dbbuild/sources/`.
As of this writing, this includes the following.

- `dbbuild/sources/alfred2018/`
- `dbbuild/sources/vandergaag2018/`
- `dbbuild/sources/staadig2019/`
- `dbbuild/sources/chen2019/`
- `dbbuild/sources/hiroaki2015/`
- `dbbuild/sources/voskoboinik2018/`
- `dbbuild/sources/puentephillips2020/`
- `dbbuild/sources/1kgp/`

Re-compiling the database to include additional data sources is a matter of creating a new directory<sup>1</sup> in `dbbuild/sources/`, populating it with appropriately formatted files (described in section "**Data Required for Each Source**"), and executing the database build workflow (the `snakemake` command described above).

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

If a particular data source provides only population data, or only marker data, etc., the other required files should simply contain an empty file except for the header line.

### `marker.tsv`

The `marker.tsv` file contains the microhaplotype marker definition.
This includes the following fields.

- `Name`: a symbol following the `mh<chromosome><initials>-<identifier>` nomenclature proposed in [Kidd (2016)](https://dx.doi.org/10.1186/s40246-016-0078-y)
- `Xref`: an identifier (or comma-separated list of identifiers) used to cross-reference the marker in other databases; this field can be empty, but if it includes one or more identifiers, each of these identifiers should be unique to this marker across all sources
- `NumVars`: the number of SNPs (or indels) that define this marker
- `Chrom`: the chromosome on which the marker is found
- `OffsetsHg37`: offsets (from the first nucleotide on the chromosome in the hg19/GRCh37 reference assembly) of the variants that define the marker, separated by commas; note that the first nucleotide of the chromosome has an offset of `0`, the second nucleotide has an offset of `1`, and so on
- `OffsetsHg38`: same as `OffsetsHg37`, but for the `GRCh38` reference assembly
- `VarRef`: comma-separated list of variant identifiers (rsIDs) associated with this marker

For example, the first few lines of the `marker.tsv` for ALFRED looks like this.

```
Name    Xref    NumVars Chrom   OffsetsHg37     OffsetsHg38     VarRef
mh17KK-014      SI664726F       3       chr17           4497060,4497088,4497096 rs333113,rs8074965,rs11657785
mh03KK-150      SI664563E       4       chr3            131927127,131927156,131927242,131927311 rs1225051,rs1225050,rs1225049,rs1225048
mh04KK-010      SI664564F       2       chr4            1985210,1985244 rs3135123,rs495367
mh04KK-011      SI664565G       2       chr4            37857268,37857332       rs6855439,rs6531591
mh04KK-013      SI664566H       5       chr4            67578383,67578461,67578473,67578538,67578583    rs13131164,rs3775866,rs11725922,rs3775867,rs17088476
mh10KK-086      SI664588L       2       chr10           95069672,95069771       rs7909236,rs17110453
mh10KK-087      SI664589M       2       chr10           105121710,105121752     rs10884095,rs1452267
mh10KK-088      SI664590E       2       chr10           133537759,133537857     rs55897648,rs2515641
mh10KK-101      SI664591F       2       chr10           133533422,133533454     rs915907,rs915908
```

> **NOTE**: We strongly recommend and prefer that a complete list of rsIDs is provided for each marker in the `VarRef` field.
> For any particular marker, if the number of rsIDs provided matches the value in `NumVars`, the database build procedure is capable of computing GRCh37 and GRCh38 offsets automatically.
> In that case, the `OffsetsHg37` and `OffsetsHg38` fields can be left blank; if they are *not* left blank, they will only be used to double-check the automatically computed offsets.
> If one or more of the SNPs defining a marker has no rsID (i.e. if the number of rsIDs is less than `NumVars`), then the `OffsetsHg37` and `OffsetsHg38` fields must be filled in and accurate: no automatic check is possible.


### `population.tsv`

The `population.tsv` file contains a description of all the populations for which microhap allele frequency data is available.
It includes the following fields.

- `ID`: a unique identifier for this population across all sources
- `Name`: a free-text description of the population, intended to be human readable
- `Xref`: optional cross-reference

For example, the first few lines of the `population.tsv` for the 1000 Genomes Project data look like this.

```
ID      Name    Xref
CHB     Han Chinese in Beijing, China   SA004058R
JPT     Japanese in Tokyo, Japan        SA004060K
CHS     Southern Han Chinese    SA004059S
CDX     Chinese Dai in Xishuangbanna, China     SA004238R
```

**NOTE**: population names need not be unique, but `Xref` identifiers must be unique across all sources.

### `frequency.tsv`

The `frequency.tsv` file contains population allele frequency data for the microhaplotype markers.
It includes the following fields.

- `Marker`: the `Name` of the marker
- `Population`: the unique identifer of the population
- `Allele`: the allele of each variant in the microhap, separated by dashes
- `Frequency`: the frequency of the allele in the specified population (a real number between 0.0 and 1.0)

For example, the first few lines of the `frequency.tsv` for van der Gaag (2018) look like this.

```
Marker  Population      Allele  Frequency
mh06PK-24844    MHDBP-383d86606a        T,C,G,C,C,C,A,A,G,A     0.0000
mh06PK-24844    MHDBP-936bc36f79        T,C,G,C,C,C,A,A,G,A     0.0000
mh06PK-24844    MHDBP-3dab7bdd14        T,C,G,C,C,C,A,A,G,A     0.1230
mh06PK-24844    MHDBP-383d86606a        T,C,G,C,C,T,A,A,G,G     0.5660
mh06PK-24844    MHDBP-936bc36f79        T,C,G,C,C,T,A,A,G,G     0.5860
mh06PK-24844    MHDBP-3dab7bdd14        T,C,G,C,C,T,A,A,G,G     0.4250
mh06PK-24844    MHDBP-383d86606a        C,C,G,C,C,C,A,A,G,A     0.0710
mh06PK-24844    MHDBP-936bc36f79        C,C,G,C,C,C,A,A,G,A     0.0000
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

For example, the entire `indel.tsv` file for ALFRED looks like this.

```
Marker	VariantIndex	Refr	Alt
mh07KK-081	0	TA	T
mh11KK-091	0	TG	T
mh22KK-064	3	AATAATT	A
```

## PermIDs

In addition to the marker names submitted by the community, MicroHapDB mints a "PermID" for each microhap marker.
These identifiers follow the pattern `MHDBM-<hash>`, where `<hash>` is computed as described below.
The PermID is designed to be stable even if reference genome assemblies are updated (thus changing the variant coordinate system), or if different researchers identify the same marker independently and assign it different names.

The hash is computed as follows.

- calculate the offset of each variant in the marker
- retrieve the sequence of the marker extent
- concatenate all of these values into a single character string, separated by commas
- compute the 32-bit [Pearson hash](pearhash/) of this string

So for example, for marker `mh07PK-38311`:

- the string is `0,6,12,58,ACCCAGAAGATACTAAGGTAAGGAAGGAGGAATTGGACTTTACTCAGAAAAGACCTCCT`
- the hash is `PearsonHash32(0,6,12,58,ACCCAGAAGATACTAAGGTAAGGAAGGAGGAATTGGACTTTACTCAGAAAAGACCTCCT) = 3ae6dc1b`

...and therefore the PermID is `MHDBM-3ae6dc1b`.

MicroHapDB uses the [pearhash](https://github.com/ze-phyr-us/pearhash) library under the MIT license to compute Pearson hashes.
