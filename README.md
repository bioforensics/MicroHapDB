# MicroHapDB

Daniel Standage, 2018  
https://github.com/bioforensics/microhapdb

**MicroHapDB** is a package designed for scientists and researchers interested in microhaplotype analysis.
This package is a distribution and convenience mechanism and does not implement any analytics itself.
All microhap data was obtained from the [Allele Frequency Database (ALFRED)][alfred]<sup>[1]</sup> at the Yale University School of Medicine.
We are not affiliated with the provisioners of this data, but we *are* enthusiatic about the utility of this data type in a variety of bioforensic applications!

## Installation

To install:

```
pip3 install microhapdb
```

To make sure the package installed correctly:

```
pip3 install pytest
py.test --pyargs microhapdb
```

MicroHapDB requires Python version 3.

## Usage

MicroHapDB provides several convenient methods to access microhaplotype data.

- a command-line interface
- a Python API
- a collection of tab-delimited text files

### Command-line interface

The following commands are available on the command line.

- `microhapdb locus`
- `microhapdb variant`
- `microhapdb allele`
- `microhapdb population`

Invoke `microhapdb <cmd> -h` to see usage instructions and examples for each command.

### Python API

Programmatic access to microhap data within Python is as simple as invoking `import microhapdb` and querying the following tables.

- `microhapdb.allelefreqs`
- `microhapdb.loci`
- `microhapdb.populations`
- `microhapdb.variants`

Each is a [Pandas][]<sup>[2]</sup> dataframe object, supporting convenient and efficient listing, subsetting, and query capabilities.
The following example demonstrates how data across the different tables can be cross-referenced.

```python
>>> import microhapdb
>>> microhapdb.loci.query('ID == "SI664558I"')
            ID        Name  Chrom      Start        End                        Variants
101  SI664558I  mh02KK-136      2  227227672  227227743  rs6714835,rs6756898,rs12617010
>>> microhapdb.populations.query('Name.str.contains("Amer")')
           ID                Name  NumChrom
41  SA000101C   African Americans       168
59  SA004047P   African Americans       122
83  SA004250L  European Americans       198
>>> microhapdb.allelefreqs.query('Locus == "SI664558I" and Population.isin(["SA000101C", "SA004047P", "SA004250L"])')
           Locus Population Allele  Frequency
38320  SI664558I  SA000101C  G,T,C      0.172
38321  SI664558I  SA000101C  G,T,A      0.103
38322  SI664558I  SA000101C  G,C,C      0.029
38323  SI664558I  SA000101C  G,C,A      0.000
38324  SI664558I  SA000101C  T,T,C      0.293
38325  SI664558I  SA000101C  T,T,A      0.063
38326  SI664558I  SA000101C  T,C,C      0.132
38327  SI664558I  SA000101C  T,C,A      0.207
38344  SI664558I  SA004047P  G,T,C      0.156
38345  SI664558I  SA004047P  G,T,A      0.148
38346  SI664558I  SA004047P  G,C,C      0.016
38347  SI664558I  SA004047P  G,C,A      0.000
38348  SI664558I  SA004047P  T,T,C      0.336
38349  SI664558I  SA004047P  T,T,A      0.049
38350  SI664558I  SA004047P  T,C,C      0.156
38351  SI664558I  SA004047P  T,C,A      0.139
38488  SI664558I  SA004250L  G,T,C      0.384
38489  SI664558I  SA004250L  G,T,A      0.202
38490  SI664558I  SA004250L  G,C,C      0.000
38491  SI664558I  SA004250L  G,C,A      0.000
38492  SI664558I  SA004250L  T,T,C      0.197
38493  SI664558I  SA004250L  T,T,A      0.000
38494  SI664558I  SA004250L  T,C,C      0.071
38495  SI664558I  SA004250L  T,C,A      0.146
```

See the [Pandas][] documentation for more details on dataframe access and query methods.

### Tab-delimited text files

The data behind MicroHapDB is contained in 4 tab-delimited text files.
If you'd prefer not to use MicroHapDB's command-line interface or Python API, it should be trivial load these files directly into R, Julia, or the data science environment of your choice.
Invoke `microhapdb files` on the command line to see the location of the installed `.tsv` files.

- `locus.tsv`: microhaplotype loci
- `variant.tsv`: variants associated with each microhap locus
- `allele.tsv`: allele frequencies for 148 loci across 84 populations
- `population.tsv`: summary of the populations studied


## Citation

If you use this package, please cite our work.

> **Standage DS** (2018) MicroHapDB: programmatic access to published microhaplotype data. GitHub repository, https://github.com/bioforensics/microhapdb.

----------

<sup>[1]</sup>Rajeevan H, Soundararajan U, Kidd JR, Pakstis AJ, Kidd KK (2012) ALFRED: an allele frequency resource for research and teaching. *Nucleic Acids Research*, 40(D1): D1010-D1015. doi:10.1093/nar/gkr924.

<sup>[2]</sup>McKinney W (2010) Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference, 51-56*.

[alfred]: https://alfred.med.yale.edu/alfred/alfredDataDownload.asp
[Pandas]: https://pandas.pydata.org
