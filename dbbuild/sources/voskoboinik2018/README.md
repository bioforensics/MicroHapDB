# Sequencing highly polymorphic haplotypes with ONT

## Citations

Voskoboinik L, Motro U, Darvasi A (2018) Facilitating complex DNA mixture interpretation by sequencing highly polymorphic haplotypes. *FSI: Genetics*, 35:136-140, [doi:10.1016/j.fsigen.2018.05.001](https://doi.org/10.1016/j.fsigen.2018.05.001).

## Build Process

The following software is required to compile the published data into the table format required by MicroHapDB.

- Python 3
- wget
- tabix (part of the [samtools/HTSlib project](https://github.com/samtools/htslib))
- [rsidx][]

A recent version of the dbSNP database for GRCh38 is also required, both the VCF file and the corresponding tabix index.

Run the following commands from the `dbbuild/sources/voskoboinik/` directory to compile the data into the table format required by MicroHapDB.

```bash
# Download 1000 Genomes Project Phase 3 data
# Note: this script is also used to build the 1kgp source; this data only needs to be downloaded once
./download.sh

# Build rsidx index if it doesn't yet exist; requires > 1 hour
rsidx index /path/to/dbSNP.vcf.gz /path/to/dbSNP.rsidx

# Compile the data tables
./compile_marker_definitions.py table1-subset.tsv /path/to/dbSNP.vcf.gz /path/to/dbSNP.rsidx
```

## Manual Pre-processing

The file `table1-subset.tsv` was created manually from Table 1 of the manuscript.


## Known Issues

Only a summary of the marker definitions is reported in the paper.
Details about each marker are not provided, and despite extensive communication with the corresponding author I was unable to resolve discrepancies between the published summary and my replication of the marker selection.

Also, the rsID `rs113012024` was merged into `rs10987426` on July 1, 2015.
The latter rsID is stored in MicroHapDB, but the former may be needed when, e.g., querying 1000 Genomes Project data.


[rsidx]: https://github.com/bioforensics/rsidx
