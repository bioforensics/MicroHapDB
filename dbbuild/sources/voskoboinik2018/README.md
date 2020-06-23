# Sequencing highly polymorphic haplotypes with ONT

## Citations

Voskoboinik L, Motro U, Darvasi A (2018) Facilitating complex DNA mixture interpretation by sequencing highly polymorphic haplotypes. *FSI: Genetics*, 35:136-140, [doi:10.1016/j.fsigen.2018.05.001](https://doi.org/10.1016/j.fsigen.2018.05.001).

## Build Process

Run the following command from the `dbbuild/sources/voskoboinik/` directory to compile the data into the table format required by MicroHapDB.

```bash
./compile_marker_definitions.py --configfile ../../config.json table1-subset.tsv
```

## Appendix

### Manual Pre-processing

The file `table1-subset.tsv` was created manually from Table 1 of the manuscript.

### Known Issues

Only a summary of the marker definitions is reported in the paper.
Details about each marker are not provided, and despite extensive communication with the corresponding author I was unable to resolve discrepancies between the published summary and my replication of the marker selection.

Also, the rsID `rs113012024` was merged into `rs10987426` on July 1, 2015.
The latter rsID is stored in MicroHapDB, but the former may be needed when, e.g., querying 1000 Genomes Project data.
