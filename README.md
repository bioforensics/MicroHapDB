[![MicroHapDB build status][travisbadge]](https://travis-ci.org/bioforensics/MicroHapDB)
[![Install with bioconda][condabadge]](http://bioconda.github.io/recipes/microhapdb/README.html)
[![BSD licensed][licensebadge]](https://github.com/bioforensics/MicroHapDB/blob/master/LICENSE.txt)

# MicroHapDB

Daniel Standage, 2018  
https://github.com/bioforensics/microhapdb

**MicroHapDB** is a package designed for scientists and researchers interested in microhaplotype analysis.
This package is a distribution and convenience mechanism and does not implement any analytics itself.
All microhaplotype data was obtained from the [Allele Frequency Database (ALFRED)][alfred]<sup>[1]</sup> at the Yale University School of Medicine and the [Leiden Open Variation Database][lovd]<sup>[2]</sup>.
We're also happy to consider integrating additional microhap data from other sources.

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
- `microhapdb.loci`
- `microhapdb.populations`
- `microhapdb.variants`

Each is a [Pandas][]<sup>[3]</sup> dataframe object, supporting convenient and efficient listing, subsetting, and query capabilities.
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
>>> import microhapdb                                                                                                  
>>> microhapdb.id_xref('mh02KK-136')                                                                                   
              ID Reference Chrom      Start        End   AvgAe  Source                                                 
128  MHDBL000129    GRCh38  chr2  227227673  227227743  3.7742  ALFRED                                                 
>>> pops = microhapdb.populations.query('Name.str.contains("Amer")')                                                   
>>> pops
             ID                Name  Source
2   MHDBP000003   African Americans  ALFRED
3   MHDBP000004   African Americans  ALFRED
21  MHDBP000022  European Americans  ALFRED
>>> f = microhapdb.frequencies                                                                                         
>>> f[(f.Locus == "MHDBL000129") & (f.Population.isin(pops.ID))]                                                       
             Locus   Population Allele  Frequency
75117  MHDBL000129  MHDBP000003  G,T,C      0.172
75118  MHDBL000129  MHDBP000003  G,T,A      0.103
75119  MHDBL000129  MHDBP000003  G,C,C      0.029
75120  MHDBL000129  MHDBP000003  G,C,A      0.000
75121  MHDBL000129  MHDBP000003  T,T,C      0.293
75122  MHDBL000129  MHDBP000003  T,T,A      0.063
75123  MHDBL000129  MHDBP000003  T,C,C      0.132
75124  MHDBL000129  MHDBP000003  T,C,A      0.207
75333  MHDBL000129  MHDBP000004  G,T,C      0.156
75334  MHDBL000129  MHDBP000004  G,T,A      0.148
75335  MHDBL000129  MHDBP000004  G,C,C      0.016
75336  MHDBL000129  MHDBP000004  G,C,A      0.000
75337  MHDBL000129  MHDBP000004  T,T,C      0.336
75338  MHDBL000129  MHDBP000004  T,T,A      0.049
75339  MHDBL000129  MHDBP000004  T,C,C      0.156
75340  MHDBL000129  MHDBP000004  T,C,A      0.139
75525  MHDBL000129  MHDBP000022  G,T,C      0.384
75526  MHDBL000129  MHDBP000022  G,T,A      0.202
75527  MHDBL000129  MHDBP000022  G,C,C      0.000
75528  MHDBL000129  MHDBP000022  G,C,A      0.000
75529  MHDBL000129  MHDBP000022  T,T,C      0.197
75530  MHDBL000129  MHDBP000022  T,T,A      0.000
75531  MHDBL000129  MHDBP000022  T,C,C      0.071
75532  MHDBL000129  MHDBP000022  T,C,A      0.146
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

<sup>[1]</sup>Rajeevan H, Soundararajan U, Kidd JR, Pakstis AJ, Kidd KK (2012) ALFRED: an allele frequency resource for research and teaching. *Nucleic Acids Research*, 40(D1): D1010-D1015, [doi:10.1093/nar/gkr924](https://doi.org/10.1093/nar/gkr924).

<sup>[2]</sup>Fokkema IF, Taschner PE, Schaafsma GC, Celli J, Laros JF, den Dunnen JT (2011) LOVD v.2.0: the next generation in gene variant databases. *Human Mutation*, 32(5): 557-63, [doi:10.1002/humu.21438](https://doi.org/10.1002/humu.21438).

<sup>[3]</sup>McKinney W (2010) Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference, 51-56*.

[alfred]: https://alfred.med.yale.edu/alfred/alfredDataDownload.asp
[lovd]: http://www.lovd.nl/3.0/home
[Pandas]: https://pandas.pydata.org
[travisbadge]: https://img.shields.io/travis/bioforensics/MicroHapDB.svg
[pypibadge]: https://img.shields.io/pypi/v/microhapdb.svg
[condabadge]: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg
[licensebadge]: https://img.shields.io/badge/license-BSD-blue.svg
