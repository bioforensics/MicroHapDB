# Panel design algorithm, mark II

## Running the procedure

- Make sure software prerequisites are installed (see below)
- Make sure databases are in place (see below)
- Edit the JSON config file `config-example.json` to point to the correct databases (or fiddle with the parameters if you're brave)
- Run the Snakemake workflow

```
snakemake --configfiles config-example.json -c 1 -p
```

### Software prerequisites

- intervaltree
- matplotlib
- networkx
- pandas
- polars
- scipy
- snakemake
- tagore
- tqdm
- upsetplot

### Required databases

- RepeatMasker track from UCSC genome browser
- dbSNP combined VCF file (along with .tbi index)

Both databases must use GRC838 coordinates.


## How it works

The design algorithm has three stages: a filtering stage, a panel scaffolding stage, and a panel fill-out stage.

1. The filtering stage applies five filters to exclude microhaps likely to perform poorly in a multiplex targeted amplicon sequencing assay: details are shown below.
2. In the scaffolding stage, each chromosome is considered independently. A *linkage graph* is constructed for each chromosome where each node represents a microhap marker and each pair of nodes is connected if the corresponding microhaps are separated by at least 9.5 Mbp (as a proxy for linkage equilibrium). All *maximal cliques* in this linkage graph are enumerated, each representing a set of mutually independent microhaps. The maximal clique with the highest aggregate Ae score is retained as the panel scaffolding for this chromosome.
3. The fill-out stage also proceeds on a chromosome-by-chromosome basis. A greedy algorithm is used to add additional independently inherited microhaps: the highest-ranked microhap by Ae that is separated by at least 9.5 Mbp from all microhaps already included in the panel is added to the panel. This is repeated until no more microhaps can be added.

### Filtering stage

The filters/masks applied during this first stage are as follows.

1. Exclude any marker within 9.5 Mbp of a forensic STR
2. Exclude any marker that overlaps with a highly conserved genomic repeat element (SINE, LINE, or LTR)
3. Exclude any marker with low-complexity sequence close to an allele-defining SNP (ADS)
4. Exclude any marker with an indel polymorphism close to an ADS
5. Exclude any marker longer than 260 bp in length

The results of each individual filter and of all aggregated filters are in `data/intermediate/`. A plot showing the number of microhaps excluded by each filter is shown in `data/results/masking-results-plot.png`.

### Whitelist

A whitelist was constructed of microhaps to include regardless of filtering status. This list primarily contained microhaps whose performance has already been demonstrated empirically in previous studies. This includes:

- All loci from ThermoFisher 74-plex
- All loci from Ken Kidd 2022 24-plex
- All loci from USC panel

It also includes two "keepers" from manual analysis of high Ae microhaps that were filtered in the preliminary stages of algorithm development (mh06SCUZJ-0528857, mh15SCUZJ-0082880).

```python
>>> table2= pd.read_csv("markers-failed-filter.tsv", sep="\t")
>>> subtable = table[(~table.FailMode.str.contains("length")) & (~table.FailMode.str.contains("str")) & (table.Ae > 9.0)].sort_values("Ae", ascending=False)
>>> subtable.to_csv("input/high-ae-filtered.tsv", sep="\t", index=False)
```


## Whence the parameters?

Some of the parameter values (declared in `config-example.json`) were selected based on informed intuition, and some based on empirical observations. The reasoning behind parameter selection is elaborated upon here.

### max_length

A limit of 300 bp is commonly used in the literature for defining a microhap. However, the Microhap Working Group elected to restrict the max extent of any microhap in a core panel to 250 bp, allowing the entire amplicon—primers and all—to fit within roughly 300 bp. So initially the value of this parameter was set to 250 bp. But during the preliminary stages of the panel design algorithm development, a handful of promising microhaps was observed right on the fence. So this parameter was marginally relaxed to 260 bp to capture a few high-value targets.

### sine/line/ltr

RepeatMasker reports a score for each genomic repeat it annotates. This score captures the extent to which any given repetitive sequence is conserved throughout the genome. Higher scores correspond to longer sequences conserved at higher fidelity, while lower scores indicate shorter sequences with more distant similarities. The distribution of these scores for all SINEs, LINEs, and LTRs in the human genome was examined, and representatives in various score ranges were observed to assess the extent of conservation (or conversely, the amount of unique sequence) corresponding to different score ranges. These observations were then used to inform the selection of cutoffs: repeat elements with scores exceeding the cutoffs were retained to filter out microhaps; repeats with lower scores were ignored at the filtering stage.

### ld_dist

A distance of 10 Mbp was initially selected as a proxy for linkage equilibrium, i.e., the physical distance required between a pair of loci to be considered independently inherited and thus suitable for the probability product rule. This threshold is applied both to filtering based on forensic STRs and to populating the linkage graph. During the preliminary stages of the panel design algorithm development, a handful of promising microhaps was observed just below that threshold. After confirming that none of these microhaps was in linkage disequilibrium with the closest candidate markers in the panel, the threshold was relaxed to 9.5 Mbp.

### max_short_mh_per_chrom

During the preliminary stages of algorithm development, it was noted that populating the linkage graph with all candidates from the chromosome would lead to maximal cliques composed of numerous short microhaps with mediocre Ae scores, crowding out microhaps with higher individual Ae scores. Rather than developing a more sophisticated clique ranking score that gives more weight to higher Ae values, it was decided to include only a handful of short microhaps in the initial linkage graph construction, and then fill in later with a greedy algorithm. This parameter limited the linkage graph to 6 microhaps < 100 bp in length for each chromosome.
