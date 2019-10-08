# MicroHapDB

[![MicroHapDB build status][travisbadge]](https://travis-ci.org/bioforensics/MicroHapDB)
[![Install with bioconda][condabadge]](http://bioconda.github.io/recipes/microhapdb/README.html)
[![BSD licensed][licensebadge]](https://github.com/bioforensics/MicroHapDB/blob/master/LICENSE.txt)

Daniel Standage, 2018-2019  
https://github.com/bioforensics/microhapdb

**MicroHapDB** is a package designed for scientists and researchers interested in microhaplotype analysis.
This package is a distribution and convenience mechanism and does not implement any analytics itself.
All microhaplotype data was obtained from public sources, including the [Allele Frequency Database (ALFRED)][alfred]<sup>[1-3]</sup> and published papers and posters<sup>[4-5]</sup>.
We're eager to consider integrating additional microhap marker and frequency data from other sources.

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

MicroHapDB requires Python version 3 and the [Pandas][] library.

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
- `microhapdb.markers`
- `microhapdb.populations`
- `microhapdb.variants`

Each is a [Pandas][]<sup>[6]</sup> dataframe object, supporting convenient and efficient listing, subsetting, and query capabilities.
There are also two auxiliary tables: one that contains a mapping of all variants to their corresponding microhap markers, and another table cross-referencing external IDs/labels/names with internal MicroHapDB identifiers.

- `microhapdb.variantmap`
- `microhapdb.idmap`

The helper function `microhapdb.id_xref` is also useful for retrieving data using any valid identifiers.
The following example demonstrates how data across the different tables can be cross-referenced.

```python
>>> import microhapdb
>>> microhapdb.id_xref('mh02KK-136')
              ID Reference Chrom      Start        End   AvgAe  Source
151  MHDBM000152    GRCh38  chr2  227227672  227227743  3.7742  ALFRED
>>> pops = microhapdb.populations.query('Name.str.contains("Amer")')
>>> pops
             ID                Name  Source
2   MHDBP000003   African Americans  ALFRED
3   MHDBP000004   African Americans  ALFRED
21  MHDBP000022  European Americans  ALFRED
>>> f = microhapdb.frequencies
>>> f[(f.Marker == "MHDBM000152") & (f.Population.isin(pops.ID))]
            Marker   Population Allele  Frequency
75117  MHDBM000152  MHDBP000003  G,T,C      0.172
75118  MHDBM000152  MHDBP000003  G,T,A      0.103
75119  MHDBM000152  MHDBP000003  G,C,C      0.029
75120  MHDBM000152  MHDBP000003  G,C,A      0.000
75121  MHDBM000152  MHDBP000003  T,T,C      0.293
75122  MHDBM000152  MHDBP000003  T,T,A      0.063
75123  MHDBM000152  MHDBP000003  T,C,C      0.132
75124  MHDBM000152  MHDBP000003  T,C,A      0.207
75333  MHDBM000152  MHDBP000004  G,T,C      0.156
75334  MHDBM000152  MHDBP000004  G,T,A      0.148
75335  MHDBM000152  MHDBP000004  G,C,C      0.016
75336  MHDBM000152  MHDBP000004  G,C,A      0.000
75337  MHDBM000152  MHDBP000004  T,T,C      0.336
75338  MHDBM000152  MHDBP000004  T,T,A      0.049
75339  MHDBM000152  MHDBP000004  T,C,C      0.156
75340  MHDBM000152  MHDBP000004  T,C,A      0.139
75525  MHDBM000152  MHDBP000022  G,T,C      0.384
75526  MHDBM000152  MHDBP000022  G,T,A      0.202
75527  MHDBM000152  MHDBP000022  G,C,C      0.000
75528  MHDBM000152  MHDBP000022  G,C,A      0.000
75529  MHDBM000152  MHDBP000022  T,T,C      0.197
75530  MHDBM000152  MHDBP000022  T,T,A      0.000
75531  MHDBM000152  MHDBP000022  T,C,C      0.071
75532  MHDBM000152  MHDBP000022  T,C,A      0.146
```

See the [Pandas][] documentation for more details on dataframe access and query methods.

### Tab-delimited text files

The data behind MicroHapDB is contained in 6 tab-delimited text files.
If you'd prefer not to use MicroHapDB's command-line interface or Python API, it should be trivial load these files directly into R, Julia, or the data science environment of your choice.
Invoke `microhapdb files` on the command line to see the location of the installed `.tsv` files.

- `marker.tsv`: microhaplotype marker definitions
- `variant.tsv`: variants associated with each microhap marker
- `allele.tsv`: allele frequencies for 148 markers across 84 populations
- `population.tsv`: summary of the populations studied
- `variantmap.tsv`: identifies which variants are associated with which markers
- `idmap.tsv`: mapping of all IDs/names/labels to internal MicroHapDB IDs


## Citation

If you use this package, please cite our work.

> **Standage DS** (2018) MicroHapDB: programmatic access to published microhaplotype data. GitHub repository, https://github.com/bioforensics/microhapdb.

----------


## References

### Published Marker collections and Allele Frequency Data

<sup>[1]</sup>Rajeevan H, Soundararajan U, Kidd JR, Pakstis AJ, Kidd KK (2012) ALFRED: an allele frequency resource for research and teaching. *Nucleic Acids Research*, 40(D1): D1010-D1015, [doi:10.1093/nar/gkr924](https://doi.org/10.1093/nar/gkr924).

<sup>[2]</sup>Kidd KK, Pakstis AJ, Speed WC, Lagace R, Wootton S, Chang J (2018) Selecting microhaplotypes optimized for different purposes. *Electrophoresis*, [doi:10.1002/elps.201800092](https://doi.org/10.1002/elps.201800092).

<sup>[3]</sup>Kidd KK, Rajeevan H (2018) ALFRED data download. *The Allele Frequency Database*, https://alfred.med.yale.edu/alfred/selectDownload/Microhap_alleleF_198.txt. Accessed December 10, 2018.

<sup>[4]</sup>van der Gaag KJ, de Leeuw RH, Laros JFJ, den Dunnen JT, de Knijff P (2018) Short hypervariable microhaplotypes: A novel set of very short high discriminating power loci without stutter artefacts. *Forensic Science International: Genetics*, 35:169-175, [doi:10.1016/j.fsigen.2018.05.008](https://doi.org/10.1016/j.fsigen.2018.05.008).

<sup>[5]</sup>Staadig A, Tillmar A (2019) Evaluation of microhaplotypes—A promising new type of forensic marker. *The 28th Congress of the International Society for Forensic Genetics*, P597.

### Supporting Software

<sup>[6]</sup>McKinney W (2010) Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference, 51-56*.


[alfred]: https://alfred.med.yale.edu/alfred/alfredDataDownload.asp
[lovd]: http://www.lovd.nl/3.0/home
[Pandas]: https://pandas.pydata.org
[travisbadge]: https://img.shields.io/travis/bioforensics/MicroHapDB.svg
[pypibadge]: https://img.shields.io/pypi/v/microhapdb.svg
[condabadge]: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg
[licensebadge]: https://img.shields.io/badge/license-BSD-blue.svg
