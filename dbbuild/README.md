# MicroHapDB Build

This directory describes the procedure used to obtain microhaplotype data from primary sources.
Most users of MicroHapDB should not have to worry about these details.
Rather, the procedure is shared in the spirit of transparency and openness necessary for independent evaluation and constructive critique.
Those who are especially curious or who wish to reproduce this work should find everything they need here.


## Prerequisites

Building MicroHapDB from scratch requires the following software and data.

- Python 3 (tested with 3.6, but will probably work with earlier 3.x versions)
- Snakemake (tested with 5.1.5)
- dbSNP database (VCF format, gzip compressed)

For convenience, let's assume for the rest of this manual that the path to the dbSNP database is stored in the environmental variable `dbsnp`.

```
export dbsnp=/path/to/the/file/dbSNP_GRCh38.vcf.gz
```

The build procedure has been tested on the Linux and Mac OS X operating systems, but OS-specific commands were intentionally avoided so the workflow should have no problems running on Windows or your cousin's obscure flavor of UNIX.


## Build step 1: data download

The first step of the procedure is to download data from the ALFRED database.
The two commands in this step are the only commands that require an Internet connection.
All subsequent commands will work offline.

```
snakemake listloci --config dbsnp=$dbsnp
snakemake fetchallloci --config dbsnp=$dbsnp
```

**Note**: Snakemake will choke on all other commands/targets if `listloci` has not yet been executed. Invoking `snakemake listloci fetchallloci ...` will not work.

**Another note**: by default, SnakeMake will only run one process at a time.
Please do not execute these initial commands in parallel mode (with the `--jobs` parameter) as that could overload the ALFRED server with numerous network requests in a very short period of time.

**Yet another note**: the `downloads/` directory is distributed with the MicroHapDB codebase as downloaded from ALFRED in September 2018.


## Build step 2: data processing

The second step of the procedure is to scrape, clean, and format the data downloaded from ALFRED, and to cross-reference the data with dbSNP.

If you have multiple processors available, you can run many steps in parallel.
For example, if you have 8 processors available add `--jobs 8` to the command below.
However, iterating over the dbSNP VCF file takes orders of magnitude more time than any other step and cannot be sped up by parallel processing, so it's a moot point.

```
snakemake tables --config dbsnp=$dbsnp
```
