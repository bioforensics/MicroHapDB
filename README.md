# MicroHapDB

[![MicroHapDB build status][travisbadge]](https://travis-ci.org/bioforensics/MicroHapDB)
[![Install with bioconda][condabadge]](http://bioconda.github.io/recipes/microhapdb/README.html)
[![BSD licensed][licensebadge]](https://github.com/bioforensics/MicroHapDB/blob/master/LICENSE.txt)

Daniel Standage, 2018-2019  
https://github.com/bioforensics/microhapdb

**MicroHapDB** is a portable database designed for scientists and researchers interested in using microhaplotype markers in forensic analysis.
The database is distributed as a collection of tabular data files in plain text, which can be queried directly or using MicroHapDB's Python API or command-line interface.
All microhaplotype marker and allele frequency data was obtained from public sources, including the [Allele Frequency Database (ALFRED)][alfred]<sup>[1-3]</sup> and published papers and posters<sup>[4-7]</sup>.
Instructions for extending your own local copy of the database with private data are provided.
However, we are eager to integrate microhap marker and frequency data from additional sources into the public database.


## Installation

For best results, install from [bioconda](https://bioconda.github.io/).

```
conda install -c bioconda microhapdb
```

To make sure the package installed correctly:

```
conda install pytest
pytest --pyargs microhapdb --doctest-modules
```

Conda ensures the correct installation of Python version 3 and the [Pandas][] library, which are required by MicroHapDB.


## Usage

MicroHapDB provides several convenient methods to query microhaplotype marker, population, and allele frequency data.

- a command-line interface
- a Python API
- a collection of tab-delimited text files

### Command-line interface

The `microhapdb` command supports four primary operations.

- `microhapdb marker`: query marker definitions using microhap names or identifiers, or variant rsIDs (such as `mh01KK-172`, `mh06PK-24844`, or `rs8192488`)
- `microhapdb population`: query population info using names or identifiers (such as `Italians` or `SA004240K`)
- `microhapdb frequency`: query allele frequencies using marker and population identifiers
- `microhapdb lookup`: query all data using any name or identifier

Invoke `microhapdb marker --help` (and so on) for query instructions and usage examples.

<img alt="MicroHapDB UNIX CLI" src="img/microhapdb-unix-cli.gif" width="600px" />

### Python API

Programmatic access to microhap data within Python is as simple as invoking `import microhapdb` and querying the following tables.

- `microhapdb.markers`
- `microhapdb.populations`
- `microhapdb.frequencies`

Each is a [Pandas][]<sup>[8]</sup> dataframe object, supporting convenient and efficient listing, subsetting, and query capabilities.

<img alt="MicroHapDB Python API" src="img/microhapdb-python-api.gif" width="600px" />

See the [Pandas][] documentation for more details on dataframe access and query methods.

MicroHapDB also includes 4 auxiliary tables, which may be useful in a variety of scenarios.

- `microhapdb.variantmap`: contains a mapping of dbSNP variants to their corresponding microhap markers
- `microhapdb.idmap`: cross-references external names and identifiers with internal MicroHapDB identifiers
- `microhapdb.sequences`: contains the sequence spanning and flanking each microhap locus
- `microhapdb.indels`: contains variant information for markers that include insertion/deletion variants

### Tab-delimited text files

The data in MicroHapDB is contained within 7 tab-delimited text files.
If you'd prefer not to use MicroHapDB's command-line interface or Python API, it should be straightfoward to load these files directly into R, Julia, or the data science environment of your choice.
Invoke `microhapdb --files` on the command line to see the location of the installed `.tsv` files.

- `frequency.tsv`: allele frequency data, indexed by marker and population
- `idmap.tsv`: mapping of external identifiers to microhap names
- `indels.tsv`: variant description for markers that include insertion/deletion variants
- `marker.tsv`: microhaplotype marker definitions
- `population.tsv`: population descriptions
- `sequences.tsv`: sequences spanning and flanking each microhap locus
- `variantmap.tsv`: a mapping of dbSNP variants to their corresponding microhap markers


## Adding Markers to MicroHapDB

> *I have some private microhap markers.
> Is it possible to include these in my MicroHapDB queries?*

Certainly!
See the [dbbuild](dbbuild/) directory for instructions on rebuilding the database with additional sources of data.

> *I have published (or am getting ready to publish) a new panel of microhap markers and allele frequencies.
> Could you add these to MicroHapDB?*

Certainly!
The instructions in the [dbbuild](dbbuild/) directory describe what data files are required.
We would be happy to assist getting data into the correct format if that would help—just let us know by opening a thread on [MicroHapDB's issue tracker](https://github.com/bioforensics/MicroHapDB/issues/new).


## Citation

If you use this package, please cite our work.

> **Standage DS** (2018) MicroHapDB: programmatic access to published microhaplotype data. GitHub repository, https://github.com/bioforensics/microhapdb.

----------


## References

### Published Marker collections and Allele Frequency Data

<sup>[1]</sup>Rajeevan H, Soundararajan U, Kidd JR, Pakstis AJ, Kidd KK (2012) ALFRED: an allele frequency resource for research and teaching. *Nucleic Acids Research*, 40(D1): D1010-D1015, [doi:10.1093/nar/gkr924](https://doi.org/10.1093/nar/gkr924).

<sup>[2]</sup>Kidd KK, Pakstis AJ, Speed WC, Lagace R, Wootton S, Chang J (2018) Selecting microhaplotypes optimized for different purposes. *Electrophoresis*, [doi:10.1002/elps.201800092](https://doi.org/10.1002/elps.201800092).

<sup>[3]</sup>Kidd KK, Rajeevan H (2018) ALFRED data download. *The Allele Frequency Database*, https://alfred.med.yale.edu/alfred/selectDownload/Microhap_alleleF_198.txt. Accessed December 7, 2018.

<sup>[4]</sup>van der Gaag KJ, de Leeuw RH, Laros JFJ, den Dunnen JT, de Knijff P (2018) Short hypervariable microhaplotypes: A novel set of very short high discriminating power loci without stutter artefacts. *Forensic Science International: Genetics*, 35:169-175, [doi:10.1016/j.fsigen.2018.05.008](https://doi.org/10.1016/j.fsigen.2018.05.008).

<sup>[5]</sup>Staadig A, Tillmar A (2019) Evaluation of microhaplotypes—A promising new type of forensic marker. *The 28th Congress of the International Society for Forensic Genetics*, P597.

<sup>[6]</sup>Hiroaki N, Fujii K, Kitayama T, Sekiguchi K, Nakanishi H, Saito K (2015) Approaches for identifying multiple-SNP haplotype blocks for use in human identification. *Legal Medicine*, 17(5):415-420, [doi:10.1016/j.legalmed.2015.06.003](https://doi.org/10.1016/j.legalmed.2015.06.003).

<sup>[7]</sup>Chen P, Deng C, Li Z, Pu Y, Yang J, Yu Y, Li K, Li D, Liang W, Zhang L, Chen F (2019) A microhaplotypes panel for massively parallel sequencing analysis of DNA mixtures. *FSI: Genetics*, 40:140-149, [doi:10.1016/j.fsigen.2019.02.018](https://doi.org/10.1016/j.fsigen.2019.02.018).

### Supporting Software

<sup>[8]</sup>McKinney W (2010) Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference, 51-56*.


[alfred]: https://alfred.med.yale.edu/alfred/alfredDataDownload.asp
[Pandas]: https://pandas.pydata.org
[travisbadge]: https://img.shields.io/travis/bioforensics/MicroHapDB.svg
[pypibadge]: https://img.shields.io/pypi/v/microhapdb.svg
[condabadge]: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg
[licensebadge]: https://img.shields.io/badge/license-BSD-blue.svg
