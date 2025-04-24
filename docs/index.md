# MicroHapDB

**MicroHapDB** is a comprehensive catalog of human microhaplotype variation.
The database integrates marker and population frequency data from numerous published research articles—see [citations](citations.md).
Allele frequency estimates and allelic diversity statistics ($A_e$) are computed for 26 global populations so that markers can be ranked and evaluated for various applications.

MicroHapDB is managed as a community resource requiring minimal infrastructure to maintain: the entire contents of the database are distributed in plain text with each copy of MicroHapDB.
The primary interface to MicroHapDB is through the console command line, which provides methods for browsing, searching, and filtering the database contents.
Full programmatic access is available via a Python API.
Alternatively, users can access the database tables directly using a spreadsheet program such as Microsoft Excel.


## Table of Contents

```{toctree}
:maxdepth: 2

starting
install
ranking
targets
extending
citations
```

If you use this database, please cite our work.

- Standage DS,  Mitchell RN (2020) MicroHapDB: A Portable and Extensible Database of All Published Microhaplotype Marker and Frequency Data. *Frontiers in Genetics* 11:781, [doi:10.3389/fgene.2020.00781](https://doi.org/10.3389/fgene.2020.00781).

MicroHapDB was created and is maintained by the Bioinformatics Group at the National Bioforensic Anaylsis Center (NBFAC).
Contact: daniel.standage@st.dhs.gov.


## Reference: command line interface

```{toctree}
:maxdepth: 1

cli-population
cli-marker
cli-frequency
cli-lookup
```


## Indices and tables

- {ref}`genindex`
- {ref}`modindex`
- {ref}`search`
