# Getting started

MicroHapDB is designed as a community resource requiring minimal infrastructure to use and maintain.
No web interface is available for MicroHapDB.
In fact, installation of the MicroHapDB software isn't strictly required to access the database contents.
For the bold and intrepid user, all database tables can be downloaded anonymously from the [https://github.com/bioforensics/MicroHapDB/tree/master/microhapdb/data](https://github.com/bioforensics/MicroHapDB/tree/master/microhapdb/data) and opened in Excel, R, Python, or any preferred data analysis environment.

The recommended mode of accessing MicroHapDB is through a single package with data and software bundled together.
Instructions for installing this package are provided [here](install.md).
Once installed, users can access the database contents in any of the following ways.


## Command line interface

MicroHapDB provides a simple and self-documenting text interface for database query and retrieval.
Using the `microhapdb` command in the terminal, a user can specify filtering criteria to select population, marker, and frequency data, and format this data in a variety of ways.

- `microhapdb population` for retrieving information on populations for which MicroHapDB has allele frequency data
- `microhapdb marker` for retrieving marker information
- `microhapdb frequency` for retrieving microhap population frequencies
- `microhapdb lookup` for retrieving individual records of any type

Executing any of the commands above with the text `--help` (e.g. `microhapdb marker --help`) will print a usage statement to the terminal with detailed instructions for configuring and running database queries.


## Python API

For users with programming experience, the contents of MicroHapDB can be accessed programmatically from the `microhapdb` Python package.
Including `import microhapdb` in the header of a Python program will provide access to the following resources.

### Database tables

MicroHapDB is comprised of three primary database tables.
In the Python API, each is stored in memory as a [Pandas](https://pandas.pydata.org) DataFrame object.

- `microhapdb.markers`
- `microhapdb.populations`
- `microhapdb.frequencies`

Additional auxiliary tables are also provided, including the following.

- `microhapdb.variantmap`: contains a mapping of dbSNP variants to their corresponding microhap markers
- `microhapdb.idmap`: cross-references external names and identifiers with internal MicroHapDB identifiers
- `microhapdb.sequences`: contains the sequence spanning and flanking each microhap locus
- `microhapdb.indels`: contains variant information for markers that include insertion/deletion variants

### Convenience functions for data retrieval

The `microhapdb.Marker` and `microhapdb.Population` modules include functions for retrieving database records based on IDs, names, genomic position, and other attributes.
Some functions return `Marker` or `Population` objects with various helpful attributes and methods.

```python
>>> marker = microhapdb.Marker.from_id("mh02USC-2pC")
>>> print(marker)
mh02USC-2pC (chr2:79025823-79025903)
>>> for marker in microhapdb.Marker.from_ids(["mh06USC-6pB", "mh17USC-17qA", "mh02USC-2pA"]):
...   print(marker)
... 
mh02USC-2pA (chr2:10810991-10811070)
mh06USC-6pB (chr6:53836249-53836261)
mh17USC-17qA (chr17:27762200-27762288)
>>> for marker in microhapdb.Marker.from_region("chr11:25000000-50000000"):
...   print(marker)
... 
mh11ZBF-001 (chr11:27379841-27379905)
mh11PK-63643 (chr11:34415814-34415851)
mh11USC-11pB (chr11:34415816-34415837)
>>> panel = list(microhapdb.Marker.from_query("Source == '10.1007/s00414-020-02483-x'"))
>>> len(panel)
59
>>> population = microhapdb.Population.from_id("IBS")
>>> print(population)
IBS     Iberian Population in Spain     1KGP
>>> for population in microhapdb.Population.from_ids(["GBR", "FIN", "CEU"]):
...   print(population)
... 
GBR     British in England and Scotland 1KGP
FIN     Finnish in Finland      1KGP
CEU     Utah Residents (CEPH) with Northern and Western European Ancestry       1KGP
>>> for population in microhapdb.Population.from_query("Name.str.contains('Afr')"):
...   print(population)
... 
MHDBP-3dab7bdd14        Africa  10.1016/j.fsigen.2018.05.008
SA000101C       African Americans       ALFRED
ACB     African Caribbeans in Barbados  1KGP
ASW     Americans of African Ancestry in SW USA 1KGP
```

A similar set of functions return DataFrame objects with subsets of the primary database tables.

```python
>>> microhapdb.Marker.table_from_ids(["mh06USC-6pB", "mh17USC-17qA", "mh02USC-2pA"])
             Name          PermID Reference  Chrom                              Offsets      Ae      In     Fst                        Source
48    mh02USC-2pA  MHDBM-1734fe04    GRCh38   chr2  10810991,10811035,10811042,10811069  2.7695  0.3702  0.2143  10.1016/j.fsigen.2019.102213
226   mh06USC-6pB  MHDBM-7cd89ff8    GRCh38   chr6           53836249,53836252,53836260  3.1711  0.1165  0.0948  10.1016/j.fsigen.2019.102213
534  mh17USC-17qA  MHDBM-cd7a9041    GRCh38  chr17  27762200,27762204,27762238,27762287  3.5538  0.0604 -0.0283  10.1016/j.fsigen.2019.102213
>>> microhapdb.Marker.table_from_region("chr11:25000000-50000000")
             Name          PermID Reference  Chrom                                            Offsets      Ae      In     Fst                        Source
368   mh11ZBF-001  MHDBM-6a26f27d    GRCh38  chr11                                  27379841,27379901  2.4158  0.0755  0.0262        10.1002/elps.201900451
369  mh11PK-63643  MHDBM-c5ce121f    GRCh38  chr11  34415814,34415816,34415818,34415835,34415836,3...     NaN     NaN     NaN  10.1016/j.fsigen.2018.05.008
370  mh11USC-11pB  MHDBM-2408c5b4    GRCh38  chr11                34415816,34415818,34415835,34415836  3.9841  0.1404  0.1346  10.1016/j.fsigen.2019.102213
>>> panel = microhapdb.Marker.table_from_query("Source == '10.1007/s00414-020-02483-x'")
>>> panel.shape
(59, 9)
>>> microhapdb.Population.table_from_ids(["GBR", "FIN", "CEU"])
      ID                                               Name Source
12   GBR                    British in England and Scotland   1KGP
26   FIN                                 Finnish in Finland   1KGP
103  CEU  Utah Residents (CEPH) with Northern and Wester...   1KGP
>>> microhapdb.Population.table_from_query("Name.str.contains('Afr')")
                 ID                                     Name                        Source
2  MHDBP-3dab7bdd14                                   Africa  10.1016/j.fsigen.2018.05.008
3         SA000101C                        African Americans                        ALFRED
4               ACB           African Caribbeans in Barbados                          1KGP
5               ASW  Americans of African Ancestry in SW USA                          1KGP
```


## Direct file access

For users that have successfully installed MicroHapDB but do not want to access the database contents through the CLI or API, database tables can be accessed directly and loaded into Excel, R, Python, or any preferred data analysis environment.
Running `microhapdb --files` on the command line will reveal the location of these files on the local system.

> *__WARNING__: Modifying the contents of the database files may cause problems with MicroHapDB. Any user wishing to sort, filter, or otherwise manipulate the contents of the core database files should instead make copies of those files and manipulate the copies.*
