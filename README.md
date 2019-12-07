# MicroHapDB

[![MicroHapDB build status][cibadge]](https://github.com/bioforensics/MicroHapDB/actions)
[![Install with bioconda][condabadge]](http://bioconda.github.io/recipes/microhapdb/README.html)
[![BSD licensed][licensebadge]](https://github.com/bioforensics/MicroHapDB/blob/master/LICENSE.txt)

Daniel Standage, 2018-2019  
https://github.com/bioforensics/microhapdb

**MicroHapDB** is a portable database intended for scientists and researchers interested in using microhaplotype markers for forensic analysis.
The database includes a comprehensive collection of marker and allele frequency data from published sources, including the [Allele Frequency Database (ALFRED)][alfred]<sup>[1-3]</sup> as well as published papers and posters<sup>[4-7]</sup>.
The entire contents of the database are distributed with each copy of MicroHapDB, and instructions for adding private data to a local copy of the database are provided.
MicroHapDB is designed to be user-friendly for both practitioners and researchers, supporting a range of access methods from browsing to simple text queries to complex queries to full programmatic access via a Python API.


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

### Browsing

Typing `microhapdb marker` on the command line will print a complete listing of all microhap markers in MicroHapDB to your terminal window.
The commands `microhapdb population` and `microhapdb frequency` will do the same for population descriptions and allele frequencies.

> *__WARNING__: it's unlikely the entire data table will fit on your screen at once, so you may have to scroll back in your terminal to view all rows of the table.*

Alternatively, the files `marker.tsv`, `population.tsv`, and `frequency.tsv` can be opened in Excel or loaded into your statistics/datascience environment of choice.
Type `microhapdb --files` on the command line to see the location of these files.

### Database queries

The `microhapdb lookup <identifier>` command searches all data tables for relevant records with a user-provided name, identifier, or description, such as `mh06PK-24844`, `rs8192488`, or `Yoruba`.

The `microhapdb marker <identifier>` command searches the microhap markers with one or more user-provided names or identifiers.
The command also supports region-based queries (such as `chr1` or `chr12:1000000-5000000`), and can print either a tabular report or a detailed report.
Run `microhapdb marker --help` for additional details.

The `microhapdb population <identifier>` command searches the population & cohort table with one or more user-provided names or identifiers.
Run `microhapdb population --help` for additional details.

The `microhapdb frequency --marker <markerID> --population <popID> --allele <allele>` command searches the allele frequency table.
The search can be restricted using all query terms (marker, population, and allele), or broadened by dropping one or more of the query terms.
Run `microhapdb frequency --help` for additional details.

<img alt="MicroHapDB UNIX CLI" src="img/microhapdb-unix-cli.gif" width="600px" />

### Python API

To access MicroHapDB from Python, simply invoke `import microhapdb` and query the following tables.

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
[cibadge]: https://github.com/bioforensics/MicroHapDB/workflows/CI%20Build/badge.svg
[pypibadge]: https://img.shields.io/pypi/v/microhapdb.svg
[condabadge]: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg
[licensebadge]: https://img.shields.io/badge/license-BSD-blue.svg
