# NRIPS/Jutendo Study

## Citations

Hiroaki N, Fujii K, Kitayama T, Sekiguchi K, Nakanishi H, Saito K (2015) Approaches for identifying multiple-SNP haplotype blocks for use in human identification. *Legal Medicine*, 17(5):415-420, [doi:10.1016/j.legalmed.2015.06.003](https://doi.org/10.1016/j.legalmed.2015.06.003).

## Build Process

Run the following command from the `dbbuild/sources/hiroaki2015/` directory to compile the data into the table format required by MicroHapDB.

```
snakemake --configfile ../../config.json -p all
```


## Appendix

### Manual Pre-processing

The file `table1-subset.tsv` was created manually from Table 1 of the manuscript.
The file `frequency.tsv` was created manually from Table 2 of the manuscript.
In both cases, marker numbers were converted into marker names using the `mh<chrom>NH-<number>` convention.

The population ID was created by appending the output of `echo $'10.1016/j.legalmed.2015.06.003\tJapanese' | md5 | cut -c 1-10` to the prefix `MHDBP-`.

### Problematic marker

The dbSNP IDs reported for marker #8 correspond to variants that span over 100kb in both GRCh37 and GRCh38.

```
RSIDs	GRCh37	GRCh38
rs1026338	69528212	68662494
rs2278918	69425762	68560044
rs2278917	69528302	68662584
		
Min	69425762	68560044
Max	69528302	68662584
Span	102541	102541
```

This marker was therefore discarded during the data processing procedure.
