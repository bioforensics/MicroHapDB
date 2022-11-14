# Marker target sequences

MicroHapDB's default tabular view for markers includes the positions (offsets) of each SNP (or indel) defined for the marker.
When designing assays for targeted amplification or hybridization enrichment of microhaplotype loci, it's important to consider the relative length of each target sequence and the location of the SNPs within the target.
MicroHapDB allows the user to specify parameters for configuring target sequences to meet desired length constraints, and recomputes SNP positions relative to the target sequence(s) for convenience.

Parameters related to target sequence definition are described below.


## Delta

```
5' delta                                                                3' delta
<---------                                                            --------->
          *                          *                               *
AGCAAACAAAGAACAGTATGTGACAGAGACTGTATCTGGTATGCAAAACTTGAAATACTTACTATCTGCCACTTTACAGA
..........G..........................G...............................C..........
..........G..........................G...............................G..........
..........T..........................C...............................C..........
..........T..........................G...............................C..........
```

The `--delta` parameter specifies the number of nucleotides beyond the first and last SNPs to extend the target sequence.
The example above shows the marker mh09USC-9qA, whose 3 SNPs occupy 60 bp from the first to the last SNP (chr9:70747687-70747747 on GRCh38), along with haplotypes observed in 1000 Genomes Project data.

In this case `--delta=10`, so 10 flanking base pairs are included on either side of the terminal SNPs.

## Minimum length

Under construction.


## Extend mode

Under construction.


## Impacted outputs

Under construction.
