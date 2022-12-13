# Wu et al., 2021

## Citations

Wu R, Li H, Li R, Peng D, Wang N, Shen X, Sun H (2021) Identification and sequencing of 59 highly polymorphic microhaplotypes for analysis of DNA mixtures. *International Journal of Legal Medicine*, 135:1137-1149, [doi:10.1007/s00414-020-02483-x](https://doi.org/10.1007/s00414-020-02483-x).


## Pre-processing

The `tables2.txt` file contains Table S2 from the manuscript's supplementary data and was exported to plain text from the published Excel file.
The `table_to_marker.py` script converts this file into the format expected by the MicroHapDB build process, stored in `marker.csv`.
The data published in the study include two dbSNP rsIDs that have been merged into other rsIDs and are no longer valid.
The script also performs the following rsID substitutions.

- rs74812635 was changed to rs602427
- rs75324027 was changed to rs602875
