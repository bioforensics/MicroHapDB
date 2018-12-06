[![MicroHapDB build status][travisbadge]](https://travis-ci.org/bioforensics/MicroHapDB)
[![PyPI version][pypibadge]](https://pypi.python.org/pypi/microhapdb)
[![BSD licensed][licensebadge]](https://github.com/bioforensics/MicroHapDB/blob/master/LICENSE.txt)

# MicroHapDB

Daniel Standage, 2018  
https://github.com/bioforensics/microhapdb

**MicroHapDB** is a package designed for scientists and researchers interested in microhaplotype analysis.
This package is a distribution and convenience mechanism and does not implement any analytics itself.
MicroHapDB is designed to work with microhap data from any source, although currently all data was obtained from the [Allele Frequency Database (ALFRED)][alfred]<sup>[1]</sup> at the Yale University School of Medicine.

## Installation

To install:

```
pip3 install microhapdb
```

To make sure the package installed correctly:

```
pip3 install pytest
pytest --pyargs microhapdb --doctest-modules
```

MicroHapDB requires Python version 3.

## Usage

MicroHapDB provides several convenient methods to access microhaplotype data.

- a command-line interface
- a Python API
- a collection of tab-delimited text files

### Command-line interface

Invoke `microhapdb --help` for a description of the command-line configuration options and several usage examples.

### Python API

Programmatic access to microhap data within Python is as simple as invoking `import microhapdb` and querying the following tables.

- `microhapdb.frequencies`
- `microhapdb.loci`
- `microhapdb.populations`
- `microhapdb.variants`

Each is a [Pandas][]<sup>[2]</sup> dataframe object, supporting convenient and efficient listing, subsetting, and query capabilities.
There are also two auxiliary tables: one that contains a mapping of all variants to their corresponding microhap loci, and another table cross-referencing external IDs/labels/names with internal MicroHapDB identifiers.

- `microhapdb.variantmap`
- `microhapdb.idmap`

The helper function `microhapdb.id_xref` is also useful for retrieving data using any valid identifiers.
The following example demonstrates how data across the different tables can be cross-referenced.

```python
>>> import microhapdb
>>> microhapdb.id_xref('mh02KK-136')
              ID Reference Chrom      Start        End  Source
182  MHDBL000183    GRCh38  chr2  227227673  227227743  ALFRED
>>> pops = microhapdb.populations.query('Name.str.contains("Amer")')
>>> pops
             ID                Name  Source
40  MHDBP000041   African Americans  ALFRED
67  MHDBP000068   African Americans  ALFRED
91  MHDBP000092  European Americans  ALFRED
>>> f = microhapdb.frequencies
>>> f[(f.Locus == "MHDBL000183") & (f.Population.isin(pops.ID))]
             Locus   Population Allele  Frequency
75117  MHDBL000183  MHDBP000041  G,T,C      0.172
75118  MHDBL000183  MHDBP000041  G,T,A      0.103
75119  MHDBL000183  MHDBP000041  G,C,C      0.029
75120  MHDBL000183  MHDBP000041  G,C,A      0.000
75121  MHDBL000183  MHDBP000041  T,T,C      0.293
75122  MHDBL000183  MHDBP000041  T,T,A      0.063
75123  MHDBL000183  MHDBP000041  T,C,C      0.132
75124  MHDBL000183  MHDBP000041  T,C,A      0.207
75333  MHDBL000183  MHDBP000068  G,T,C      0.156
75334  MHDBL000183  MHDBP000068  G,T,A      0.148
75335  MHDBL000183  MHDBP000068  G,C,C      0.016
75336  MHDBL000183  MHDBP000068  G,C,A      0.000
75337  MHDBL000183  MHDBP000068  T,T,C      0.336
75338  MHDBL000183  MHDBP000068  T,T,A      0.049
75339  MHDBL000183  MHDBP000068  T,C,C      0.156
75340  MHDBL000183  MHDBP000068  T,C,A      0.139
75525  MHDBL000183  MHDBP000092  G,T,C      0.384
75526  MHDBL000183  MHDBP000092  G,T,A      0.202
75527  MHDBL000183  MHDBP000092  G,C,C      0.000
75528  MHDBL000183  MHDBP000092  G,C,A      0.000
75529  MHDBL000183  MHDBP000092  T,T,C      0.197
75530  MHDBL000183  MHDBP000092  T,T,A      0.000
75531  MHDBL000183  MHDBP000092  T,C,C      0.071
75532  MHDBL000183  MHDBP000092  T,C,A      0.146
```

See the [Pandas][] documentation for more details on dataframe access and query methods.

### Tab-delimited text files

The data behind MicroHapDB is contained in 6 tab-delimited text files.
If you'd prefer not to use MicroHapDB's command-line interface or Python API, it should be trivial load these files directly into R, Julia, or the data science environment of your choice.
Invoke `microhapdb files` on the command line to see the location of the installed `.tsv` files.

- `locus.tsv`: microhaplotype loci
- `variant.tsv`: variants associated with each microhap locus
- `allele.tsv`: allele frequencies for 148 loci across 84 populations
- `population.tsv`: summary of the populations studied
- `variantmap.tsv`: shows which variants are associated with which loci
- `idmap.tsv`: mapping of all IDs/names/labels to internal MicroHapDB IDs


## Citation

If you use this package, please cite our work.

> **Standage DS** (2018) MicroHapDB: programmatic access to published microhaplotype data. GitHub repository, https://github.com/bioforensics/microhapdb.

----------

<sup>[1]</sup>Rajeevan H, Soundararajan U, Kidd JR, Pakstis AJ, Kidd KK (2012) ALFRED: an allele frequency resource for research and teaching. *Nucleic Acids Research*, 40(D1): D1010-D1015. doi:10.1093/nar/gkr924.

<sup>[2]</sup>McKinney W (2010) Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference, 51-56*.

[alfred]: https://alfred.med.yale.edu/alfred/alfredDataDownload.asp
[Pandas]: https://pandas.pydata.org
[travisbadge]: https://img.shields.io/travis/bioforensics/MicroHapDB.svg
[pypibadge]: https://img.shields.io/pypi/v/microhapdb.svg
[licensebadge]: https://img.shields.io/badge/license-BSD-blue.svg
