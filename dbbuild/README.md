# MicroHapDB Build

This directory describes the procedure used to obtain microhaplotype data from primary sources.
Most users of MicroHapDB should not have to worry about these details.
Rather, the procedure is shared in the spirit of transparency and openness necessary for independent evaluation and constructive critique.
Those who are especially curious or who wish to reproduce this work should find everything they need here.


## Prerequisites

Building MicroHapDB from scratch requires the following software and data.

- Python 3 (tested with 3.6, but will probably work with earlier 3.x versions)
- Pandas (tested with 0.23.4)
- Snakemake (tested with 5.1.5)
- dbSNP database (VCF format, gzip compressed)
- tabix (tested with 1.9)

For convenience, let's assume for the rest of this manual that the path to the dbSNP database is stored in the environmental variable `dbsnp`.

```
export dbsnp=/path/to/the/file/dbSNP_GRCh38.vcf.gz
```

The build procedure has been tested on the Linux and Mac OS X operating systems.


## Data download

> **NOTE**: This step should not need to be repeated and is described only for the purpose of full disclosure.

Important variant information for ALFRED microhaplotypes are unavailable in convenient summary form.
The `ALFRED.Snakefile` workflow was used to retrieve HTML summary pages for microhap loci over the network from the ALFRED database.

```
snakemake --snakefile ALFRED.Snakefile -p loci
```

These files are stored in the `alfred/downloads/locus-detail/` directory for scraping by the main data processing workflow.
All other data required for the build has been downloaded as described in the `alfred/` and `lovd/` directories.

> **ANOTHER NOTE**: As of December 2018, funding for the ALFRED database is set to expire.
> It is uncertain how long this build step will work.
> Fortunately the data have been captured and stored here to enable future reference.


## Data processing

The main build procedure is implemented in `Snakefile` and scrapes, cleans, formats, and cross-references data from ALFRED, LOVD, and dbSNP.
It is invoked like so.

```
snakemake tables --config dbsnp=$dbsnp
```

To speed up the build with multiple processes:

```
snakemake tables --cores 8 --config dbsnp=$dbsnp
```
