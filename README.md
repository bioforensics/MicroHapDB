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
              ID Reference Chrom      Start        End   AvgAe  Source
128  MHDBL000129    GRCh38  chr2  227227672  227227743  3.7742  ALFRED
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


## References

### Variant Databases

<sup>[1]</sup>Rajeevan H, Soundararajan U, Kidd JR, Pakstis AJ, Kidd KK (2012) ALFRED: an allele frequency resource for research and teaching. *Nucleic Acids Research*, 40(D1): D1010-D1015, [doi:10.1093/nar/gkr924](https://doi.org/10.1093/nar/gkr924).

<sup>[2]</sup>Fokkema IF, Taschner PE, Schaafsma GC, Celli J, Laros JF, den Dunnen JT (2011) LOVD v.2.0: the next generation in gene variant databases. *Human Mutation*, 32(5): 557-63, [doi:10.1002/humu.21438](https://doi.org/10.1002/humu.21438).

### Supporting Software

<sup>[3]</sup>McKinney W (2010) Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference, 51-56*.

### Microhaplotype Studies

Bulbul O, Pakstis AJ, Soundararajan U, Gurkan C, Brissenden JE, Roscoe JM, Evsanaa B, Togtokh A, Paschou P, Grigorenko EL, Gurwitz D, Wootton S, Lagace R, Chang J, Speed WC, Kidd KK (2018) Ancestry inference of 96 population samples using microhaplotypes. *International Journal of Legal Medicine*, **132**:703-711, [doi:10.1007/s00414-017-1748-6](https://doi.org/10.1007/s00414-017-1748-6).

Chen P, Yin C, Li Z, Pu Y, Yu Y, Zhao P, Chen D, Liang W, Zhang L, Chen F (2018) Evaluation of the Microhaplotypes panel for DNA mixture analyses. *Forensic Science International: Genetics*, **35**:149-155, [doi:10.1016/j.fsigen.2018.05.003](https://doi.org/10.1016/j.fsigen.2018.05.003).

Chen P, Zhu W, Tong F, Pu Y, Yu Y, Huang S, Li Z, Zhang L, Liang W, Chen F (2018) Identifying novel microhaplotypes for ancestry inference. *International Journal of Legal Medicine*, [doi:10.1007/s00414-018-1881-x](https://doi.org/10.1007/s00414-018-1881-x).

Kidd KK, Pakstis AJ, Speed WC, Lagace R, Wootton S, Chang J (2018) Selecting microhaplotypes optimized for different purposes. *Electrophoresis*, [doi:10.1002/elps.201800092](https://doi.org/10.1002/elps.201800092).

Kidd KK, Speed WC, Pakstis AJ, Podini DS, Lagace R, Chang J, Wootton S, Haigh E, Soundararajan U (2017) Evaluating 130 Microhaplotypes across a Global Set of 83 Populations. *Forensic Science International: Genetics*, **29**:29-37, [doi:10.1016/j.fsigen.2017.03.014](https://doi.org/10.1016/j.fsigen.2017.03.014).

Kidd KK, Pakstis AJ, Speed WC, Lagace R, Chang J, Wootton S, Haigh E, Kidd JR (2014) Current sequencing technology makes microhaplotypes a powerful new type of genetic marker for forensics. *Forensic Science International: Genetics*, **12**:215-224, [doi:10.1016/j.fsigen.2014.06.014](https://doi.org/10.1016/j.fsigen.2014.06.014)

Kidd KK, Pakstis AJ, Speed WC, Lagace R, Chang J, Wootton S, Ihuegbu N (2013) Microhaplotype loci are a powerful new type of forensic marker". *Forensic Science International: Genetics supplement series*, **4**:e123-124, [doi:10.1016/j.fsigss.2013.10.063](https://doi.org/10.1016/j.fsigss.2013.10.063).

van der Gaag KJ, de Leeuw RH, Laros JFJ, den Dunnen JT, de Knijff P (2018) Short hypervariable microhaplotypes: A novel set of very short high discriminating power loci without stutter artefacts. *Forensic Science International: Genetics*, 35:169-175, [doi:10.1016/j.fsigen.2018.05.008](https://doi.org/10.1016/j.fsigen.2018.05.008).


[alfred]: https://alfred.med.yale.edu/alfred/alfredDataDownload.asp
[lovd]: http://www.lovd.nl/3.0/home
[Pandas]: https://pandas.pydata.org
[travisbadge]: https://img.shields.io/travis/bioforensics/MicroHapDB.svg
[pypibadge]: https://img.shields.io/pypi/v/microhapdb.svg
[condabadge]: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg
[licensebadge]: https://img.shields.io/badge/license-BSD-blue.svg
