# NYGC 1000 Genomes Project Update

## Citations

Byrska-Bishop M, Evani US, Zhao X, et al. (2022) High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios. *Cell* 185:18, [doi:10.1016/j.cell.2022.08.004](https://doi.org/10.1016/j.cell.2022.08.004).


## Build process

If `dbbuild/marker.csv` has been updated, copy that file to `marker-latest.csv` in this directory.

```
cp ../../marker.csv marker-latest.csv
```

Then haplotype frequencies and A<sub>e</sub> values can be recomputed as follows.

```
snakemake -c 8 --config dir_1kgp=databases/1000Genomes/ refr=databases/hg38.fasta -p
```


## Required databases

The data from this study to too large to bundle with the MicroHapDB git repository.
Prior to (re-)building this data set, all VCF and TBI files must be downloaded from the following URL: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/.
The files can be placed anywhere on your system, so long as the directory path is correctly specified in the `1kgp_dir` configuration as shown above.

The GRCh38 human reference genome assembly should also be downloaded from the url https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz, decompressed, and specified in the `refr` configuration as shown above.
