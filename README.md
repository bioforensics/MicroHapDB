# MicroHapDB

[![MicroHapDB build status][cibadge]](https://github.com/bioforensics/MicroHapDB/actions)
[![Install with bioconda][condabadge]](http://bioconda.github.io/recipes/microhapdb/README.html)
[![BSD licensed][licensebadge]](https://github.com/bioforensics/MicroHapDB/blob/master/LICENSE.txt)

NBFAC, 2018-2022
https://github.com/bioforensics/microhapdb

**MicroHapDB** is a portable database intended for scientists and researchers interested in microhaplotypes for forensic analysis.
The database includes a comprehensive collection of marker and allele frequency data from numerous databases and published research articles.<sup>[5-19]</sup>
Effective number of allele (*A<sub>e</sub>*)<sup>[2]</sup> and informativeness for assignment (*I<sub>n</sub>*)<sup>[3]</sup> statistics are included so that markers can be ranked for different forensic applications.
The entire contents of the database are distributed with each copy of MicroHapDB, and instructions for adding private data to a local copy of the database are provided.
MicroHapDB is designed to be user-friendly for both practitioners and researchers, supporting a range of access methods from browsing and simple text queries to complex queries and full programmatic access via a Python API.
MicroHapDB is also designed as a community resource requiring minimal infrastructure to use and maintain.


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

There is no web interface for MicroHapDB.
The database must be installed locally to the computer(s) on which it will be used.
The entire contents of the database are distributed with each copy of MicroHapDB.
Once installed, users can access the database contents in any of the following ways.
(Click on each arrow for more information.)

<details>
    <summary>Command line interface</summary>

### Command line interface

MicroHapDB provides a simple and user-friendly text interface for database query and retrieval.
Using the `microhapdb` command in the terminal, a user can provide filtering criteria to select population, marker, and frequency data, and format this data in a variety of ways.

- `microhapdb population` for retrieving information on populations for which MicroHapDB has allele frequency data
- `microhapdb marker` for retrieving marker information
- `microhapdb frequency` for retrieving microhap population frequencies
- `microhapdb lookup` for retrieving individual records of any type

Executing `microhapdb population --help`, `microhapdb marker --help`, and so on in the terminal will print a usage statement with detailed instructions for configuring and running database queries.

</details>

<details>
    <summary>Python API</summary>

### Python API

For users with programming experience, the contents of MicroHapDB can be accessed programmatically from the `microhapdb` Python package.
Including `import microhapdb` in the header of a Python program will provide access to the following resources.

- database tables, stored in memory as [Pandas][]<sup>[1]</sup> DataFrame objects
    - `microhapdb.markers`
    - `microhapdb.populations`
    - `microhapdb.frequencies`
- convenience functions for data retrieval
    - some return `Marker` or `Population` objects with numerous auxiliary attributes and methods
        - `marker = microhapdb.Marker.from_id("mh02USC-2pC")`
        - `for marker in microhapdb.Marker.from_ids(idlist):`
        - `for marker in microhapdb.Marker.from_query("Source == '10.1007/s00414-020-02483-x'"):`
        - `for marker in microhapdb.Marker.from_region("chr11:25000000-50000000"):`
        - `population = microhapdb.Population.from_id("IBS")`
        - `for population in microhapdb.Population.from_ids(idlist):`
        - `for population in microhapdb.Population.from_query("Name.str.contains('Afr')"):`
    - others return DataFrame objects with subsets of the primary database tables
        - `result = microhapdb.Marker.table_from_ids(idlist)`
        - `result = microhapdb.Marker.table_from_query("Source == '10.1007/s00414-020-02483-x'")`
        - `result = microhapdb.Marker.table_from_region("chr11:25000000-50000000")`
        - `result = microhapdb.Population.table_from_ids(idlist):`
        - `result = microhapdb.Population.table_from_query("Name.str.contains('Afr')"):`
</details>

<details>
    <summary>Direct file access</summary>

### Direct file access

For users that want to circumvent MicroHapDB's CLI and API, the database tables can be accessed directly and loaded into R, Python, Excel, or any preferred environment.
Running `microhapdb --files` on the command line will reveal the location of these files.

> *__WARNING__: Modifying the contents of the database files may cause problems with MicroHapDB. Any user wishing to sort, filter, or otherwise manipulate the contents of the core database files should instead copy those files and manipulate the copies.*
</details>


## Ranking Markers

MicroHapDB provides three criteria for ranking markers.

- `Ae`: the marker's effective number of alleles (*A<sub>e</sub>*) computed individually for 26 populations<sup>[2]</sup>; by default, the average of the 26 populations is shown, but the `--ae-pop` flag or the `microhapdb.set_ae_population` function can be used to specify a single population for with to display A<sub>e</sub> values
- `In`: Rosenberg's informativeness for assignment (*I<sub>n</sub>*) computed on 26 populations<sup>[3]</sup>
- `Fst`: fixation index (*F<sub>ST</sub>*) computed on 26 populations<sup>[4]</sup>

The A<sub>e</sub> statistic is a measure of the *within-population* allelic variation at a locus, which corresponds to the marker's diagnostic power for identification purposes.
The I<sub>n</sub> and F<sub>ST</sub> statistics measure *between-population* allelic variation at a locus, which corresponds to the marker's utility for predicting population of origin.
Phased genotypes for 2,504 individuals from Phase 3 of the 1000 Genomes Project<sup>[14]</sup> are used to calculate these statistics.


## Adding Markers to MicroHapDB

> *I have some private microhap markers.
> Is it possible to include these in my MicroHapDB queries without making them public?*

Certainly!
See the [dbbuild](dbbuild/) directory for instructions on updating a private copy the database with additional sources of data.

> *I have published (or am getting ready to publish) a new panel of microhap markers and allele frequencies.
> Could you add these to MicroHapDB?*

Certainly!
The instructions in the [dbbuild](dbbuild/) directory describe what data files are required.
We would be happy to assist getting data into the correct format if that would help—just let us know by opening a thread on [MicroHapDB's issue tracker](https://github.com/bioforensics/MicroHapDB/issues/new).

> *My favorite data set of population-level genomic variation is not available in MicroHapDB.
> Could you add it?*

We would be happy to consider including any data set with compatible licensing or user agreements.
The [dbbuild](dbbuild/) directory describes the information required for marker and population definitions.
Computing allele frequencies from population data typically requires variant calls (VCFs) with phased genotypes for individuals.
Use [MicroHapDB's issue tracker](https://github.com/bioforensics/MicroHapDB/issues/new) to contact us with any questions about including new public data sets.


## Citation

If you use this database, please cite our work.

> Standage DS,  Mitchell RN (2020) MicroHapDB: A Portable and Extensible Database of All Published Microhaplotype Marker and Frequency Data. *Frontiers in Genetics* 11:781, [doi:10.3389/fgene.2020.00781](https://doi.org/10.3389/fgene.2020.00781).

MicroHapDB was created and is maintained by the Bioinformatics Group at the National Bioforensic Anaylsis Center (NBFAC).


----------


## References

### Supporting Software

<sup>[1]</sup>McKinney W (2010) Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference, 51-56*.

### Ranking Statistics

<sup>[2]</sup>Crow JF, Kimura M (1970) <u>An Introduction to Population Genetics Theory</u>. New York, Harper & Row.

<sup>[3]</sup>Rosenberg NA, Li LM, Ward R, Pritchard JK (2003) Informativeness of genetic markers for inference of ancestry. *American Journal of Human Genetics*, 73(6):1402–1422, [doi:10.1086/380416](https://doi.org/10.1086/380416).

<sup>[4]</sup>Weir B, Cockerham, C (1984) Estimating F-Statistics for the Analysis of Population Structure. *Evolution* 38(6):1358-1370, [doi:10.2307/2408641](https://doi.org/10.2307/2408641).

### Published Marker collections and Allele Frequency Data

<sup>[5]</sup>Rajeevan H, Soundararajan U, Kidd JR, Pakstis AJ, Kidd KK (2012) ALFRED: an allele frequency resource for research and teaching. *Nucleic Acids Research*, 40(D1): D1010-D1015, [doi:10.1093/nar/gkr924](https://doi.org/10.1093/nar/gkr924).

<sup>[6]</sup>Kidd KK, Pakstis AJ, Speed WC, Lagace R, Wootton S, Chang J (2018) Selecting microhaplotypes optimized for different purposes. *Electrophoresis*, [doi:10.1002/elps.201800092](https://doi.org/10.1002/elps.201800092).

<sup>[7]</sup>Kidd KK, Rajeevan H (2018) ALFRED data download. *The Allele Frequency Database*, https://alfred.med.yale.edu/alfred/selectDownload/Microhap_alleleF_198.txt. Accessed December 7, 2018.

<sup>[8]</sup>van der Gaag KJ, de Leeuw RH, Laros JFJ, den Dunnen JT, de Knijff P (2018) Short hypervariable microhaplotypes: A novel set of very short high discriminating power loci without stutter artefacts. *Forensic Science International: Genetics*, 35:169-175, [doi:10.1016/j.fsigen.2018.05.008](https://doi.org/10.1016/j.fsigen.2018.05.008).

<sup>[9]</sup>Staadig A, Tillmar A (2019) Evaluation of microhaplotypes—A promising new type of forensic marker. *The 28th Congress of the International Society for Forensic Genetics*, P597.

<sup>[10]</sup>Hiroaki N, Fujii K, Kitayama T, Sekiguchi K, Nakanishi H, Saito K (2015) Approaches for identifying multiple-SNP haplotype blocks for use in human identification. *Legal Medicine*, 17(5):415-420, [doi:10.1016/j.legalmed.2015.06.003](https://doi.org/10.1016/j.legalmed.2015.06.003).

<sup>[11]</sup>Chen P, Deng C, Li Z, Pu Y, Yang J, Yu Y, Li K, Li D, Liang W, Zhang L, Chen F (2019) A microhaplotypes panel for massively parallel sequencing analysis of DNA mixtures. *FSI: Genetics*, 40:140-149, [doi:10.1016/j.fsigen.2019.02.018](https://doi.org/10.1016/j.fsigen.2019.02.018).

<sup>[12]</sup>Voskoboinik L, Motro U, Darvasi A (2018) Facilitating complex DNA mixture interpretation by sequencing highly polymorphic haplotypes. *FSI: Genetics*, 35:136-140, [doi:10.1016/j.fsigen.2018.05.001](https://doi.org/10.1016/j.fsigen.2018.05.001).

<sup>[13]</sup>de la Puente M, Phillips C, Xavier C, Amigo J, Carracedo A, Parson W, Lareu MV (2020) Building a custom large-scale panel of novel microhaplotypes for forensic identification using MiSeq and Ion S5 massively parallel sequencing systems. *FSI: Genetics*, 45:102213, [doi:10.1016/j.fsigen.2019.102213](https://doi.org/10.1016/j.fsigen.2019.102213).

<sup>[14]</sup>Auton A, Abecasis G, Altshuler D, et al. (2015) A global reference for human genetic variation. *Nature* 526:68–74, [doi:10.1038/nature15393](https://doi.org/10.1038/nature15393).

<sup>[15]</sup>Gandotra N, Speed WC, Qin W, Tang Y, Pakstis AJ, Kidd KK, Scharfe C (2020) Validation of novel forensic DNA markers using multiplex microhaplotype sequencing. *Forensic Science International: Genetics*, **47**:102275, [doi:10.1016/j.fsigen.2020.102275](https://doi.org/10.1016/j.fsigen.2020.102275).

<sup>[16]</sup>Sun S, Liu Y, Li J, Yang Z, Wen D, Liang W, Yan Y, Yu H, Cai J, Zha L (2020) Development and application of a nonbinary SNP-based microhaplotype panel for paternity testing involving close relatives. *FSI: Genetics*, 46:102255, [doi:10.1016/j.fsigen.2020.102255](https://doi.org/10.1016/j.fsigen.2020.102255).

<sup>[17]</sup>Kureshi A, Li J, Wen D, Sun S, Yang Z, Zha L (2020) Construction and forensic application of 20 highly polymorphic microhaplotypes. *Royal Society Open Science*, **7**(5):191937, [doi:10.1098/rsos.191937](https://doi.org/10.1098/rsos.191937).

<sup>[18]</sup>Jin XY, Cui W, Chen C, Guo YX, Zhang XR, Xing GH, Lan JW, Zhu BF (2020) Developing and population analysis of a new multiplex panel of 18 microhaplotypes and compound markers using next generation sequencing and its application in the Shaanxi Han population. *Electrophoresis*, **41**(13-14):1230-1237, [doi:10.1002/elps.201900451](https://doi.org/10.1002/elps.201900451).

<sup>[19]</sup>Wu R, Li H, Li R, Peng D, Wang N, Shen X, Sun H (2021) Identification and sequencing of 59 highly polymorphic microhaplotypes for analysis of DNA mixtures. *International Journal of Legal Medicine*, 135:1137-1149, [doi:10.1007/s00414-020-02483-x](https://doi.org/10.1007/s00414-020-02483-x).

[alfred]: https://alfred.med.yale.edu/alfred/alfredDataDownload.asp
[Pandas]: https://pandas.pydata.org
[cibadge]: https://github.com/bioforensics/MicroHapDB/workflows/CI%20Build/badge.svg
[pypibadge]: https://img.shields.io/pypi/v/microhapdb.svg
[condabadge]: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg
[licensebadge]: https://img.shields.io/badge/license-BSD-blue.svg
