# Short hypervariable microhaplotypes (van der Gaag *et al.* 2018)

## Citations

van der Gaag KJ, de Leeuw RH, Laros JFJ, den Dunnen JT, de Knijff P (2018) Short hypervariable microhaplotypes: A novel set of very short high discriminating power loci without stutter artefacts. *Forensic Science International: Genetics*, 35:169-175, [doi:10.1016/j.fsigen.2018.05.008](https://doi.org/10.1016/j.fsigen.2018.05.008).

## Build Process

Run the following command from the `dbbuild/sources/vandergaag2018/` directory to compile the data into the table format required by MicroHapDB.

```bash
snakemake --cores 1
```


## Appendix

### Manual Pre-processing

The file `figure-S1.txt` was manually created using the supplementary PDF file from the paper (`original/mmc1.pdf`).
Genomic coordinates of each marker were determined by copying the sequence from the PDF, editing in some cases, and pasting into a UCSC BLAT search.
The offset of each variant site from the first nucleotide in the marker was manually checked, and random spot checking was used to verify the accuracy.

The file `marker-rsids.tsv` was created by performing manual lookups of the variant positions and determining the correct rsID from the variant annotation.
This task was performed by two folks independently (@standage and @rnmitchell) and the results were cross-checked for accuracy.

The file `table-S2.txt` was created by copying the table from the supplementary XLSX file from the paper (`original/mmc4.xlsx`) and pasting into a text editor.

The file `population.tsv` was created manually.

The population IDs were created by appending the output of the following commands to the prefix `MHDBP-`.

```
$ echo $'10.1016/j.fsigen.2018.05.008\tAfrica' | md5 | cut -c 1-10
3dab7bdd14
$ echo $'10.1016/j.fsigen.2018.05.008\tAsia' | md5 | cut -c 1-10
936bc36f79
$ echo $'10.1016/j.fsigen.2018.05.008\tNL' | md5 | cut -c 1-10
383d86606a
```
