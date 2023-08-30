# MicroHapDB

[![MicroHapDB build status][cibadge]](https://github.com/bioforensics/MicroHapDB/actions)
[![Install with bioconda][condabadge]](http://bioconda.github.io/recipes/microhapdb/README.html)
[![BSD licensed][licensebadge]](https://github.com/bioforensics/MicroHapDB/blob/master/LICENSE.txt)

NBFAC, 2018-2023
https://github.com/bioforensics/microhapdb

**MicroHapDB** is a comprehensive database of microhaplotype data for human forensics.
The database integrates marker and population frequency data from numerous published research articles.
Polymorphism statistics ($A_e$) are computed for 26 global populations so that markers can be ranked for various forensic applications.

MicroHapDB is designed to be user-friendly for researchers, technology developers, and forensic practitioners, supporting a range of access methods from browsing and simple text queries to complex queries and full programmatic access via a Python API.
It is designed as a community resource requiring minimal infrastructure to use and maintain: the entire contents of the database are distributed in plain text with each copy of MicroHapDB.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/bioforensics/MicroHapDB/master?labpath=binder%2Fdemo_v0.9.ipynb)


## Installation

See [this page](https://microhapdb.readthedocs.io/en/latest/install.html) for complete installation instructions.
If this isn't your first rodeo, the commands below are provided as a quick reference.

```
conda install -c bioconda microhapdb pytest
pytest --pyargs microhapdb --doctest-modules
```

## Usage

MicroHapDB provides several methods to access the contents of a locally installed database.
The [MicroHapDB documentation](https://microhapdb.readthedocs.io/) includes a ["Getting started" guide](https://microhapdb.readthedocs.io/en/latest/starting.html) as well as a comprehensive reference for running MicroHapDB on the command line.

> *A particularly intrepid and curious user may also be so bold as to download the core database tables directly from [GitHub](https://github.com/bioforensics/MicroHapDB/tree/master/microhapdb/data).*

MicroHapDB includes [statistics for ranking markers](https://microhapdb.readthedocs.io/en/latest/ranking.html), tools for [panel design](targets), and instructions for [adding markers to a private local copy of the database](https://microhapdb.readthedocs.io/en/latest/extending.html).


## Citation

If you use this database, please cite our work.

> Standage DS,  Mitchell RN (2020) MicroHapDB: A Portable and Extensible Database of All Published Microhaplotype Marker and Frequency Data. *Frontiers in Genetics* 11:781, [doi:10.3389/fgene.2020.00781](https://doi.org/10.3389/fgene.2020.00781).

MicroHapDB was created and is maintained by the Bioinformatics Group at the National Bioforensic Anaylsis Center (NBFAC).

Additional references are available [on this page](https://microhapdb.readthedocs.io/en/latest/install.html).


[alfred]: https://alfred.med.yale.edu/alfred/alfredDataDownload.asp
[Pandas]: https://pandas.pydata.org
[cibadge]: https://github.com/bioforensics/MicroHapDB/workflows/CI%20Build/badge.svg
[pypibadge]: https://img.shields.io/pypi/v/microhapdb.svg
[condabadge]: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg
[licensebadge]: https://img.shields.io/badge/license-BSD-blue.svg
