# MicroHapDB

**MicroHapDB** is a portable database intended for scientists and researchers interested in microhaplotypes for forensic analysis.
The database includes a comprehensive collection of marker and allele frequency data from numerous databases and published research articles—see [citations](citations.md).
Effective number of allele (*Ae*) and informativeness for assignment (*In*) statistics are included so that markers can be ranked for different forensic applications.
The entire contents of the database are distributed with each copy of MicroHapDB, and instructions for adding private data to a local copy of the database are provided.
MicroHapDB is designed to be user-friendly for both practitioners and researchers, supporting a range of access methods from browsing and simple text queries to complex queries and full programmatic access via a Python API.
MicroHapDB is also designed as a community resource requiring minimal infrastructure to use and maintain.


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
