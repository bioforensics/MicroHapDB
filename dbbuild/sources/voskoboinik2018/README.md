# NRIPS/Jutendo Study

## Citations

Voskoboinik L, Motro U, Darvasi A (2018) Facilitating complex DNA mixture interpretation by sequencing highly polymorphic haplotypes. *FSI: Genetics*, 35:136-140, [doi:10.1016/j.fsigen.2018.05.001](https://doi.org/10.1016/j.fsigen.2018.05.001).

## Build Process

The following software is required to compile the published data into the table format required by MicroHapDB.

- Python 3
- tabix (part of the [samtools/HTSlib project](https://github.com/samtools/htslib))
- [rsidx][]

A recent version of the dbSNP database for GRCh38 is also required, both the VCF file and the corresponding tabix index.

Run the following commands from the `dbbuild/sources/voskoboinik/` directory to compile the data into the table format required by MicroHapDB.

```bash
# Build rsidx index if it doesn't yet exist; requires > 1 hour
rsidx index /path/to/dbSNP.vcf.gz /path/to/dbSNP.rsidx

# Compile the data tables
./compile_marker_definitions.py table1-subset.tsv /path/to/dbSNP.vcf.gz /path/to/dbSNP.rsidx
```

## Manual Pre-processing

The file `table1-subset.tsv` was created manually from Table 1 of the manuscript.



[rsidx]: https://github.com/bioforensics/rsidx
