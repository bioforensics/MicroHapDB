# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## Unreleased

### Added 
- A new table flagging markers with repetitive content (see #147).


## [0.10.1] 2023-10-13

### Fixed
- Bug with offsets table (`marker --format=offsets`) when multiple markers are defined for a locus (#144).


## [0.10] 2023-09-15

### Added
- Marker definitions for NimaGen 29-plex (see #134).
- 50 marker definitions from Zhang et al 2023 (see #136).
- 191 marker definitions from Du et al 2023 (see #137).
- Jupyter notebook for plotting database growth over time (see #135).

### Changed
- Cleanup of database build code (see #133).
- Improved marker detail view now includes locus, marker, synonym, and source info (see #138).
- Updated documentation (see #139, #142).

### Fixed
- Off-by-one errors in Fan 2022 marker definitions (see #140).


## [0.9] 2023-07-03

### Added
- Italian population frequencies for 98 markers from Turchi et al 2019 (see #107).
- 21 marker definitions and population frequences from Zou et al 2022 (see #109).
- Check for overlapping microhaps from different sources during the database build (#112).
- Over 1,000 marker definitions from four whole-genome screens from Yu et al 2022 (see #115, #118, #128).
- New `nomenclature` module for parsing and validating microhap identifiers (see #119).
- New `summarize` subcommand for summarizing MicroHapDB database contents (see #120).
- 1,000 marker definitions from Zhu et al 2023 (see #125, #131).
- New check during database build time to make sure that markers sharing the same name actually occupy the same locus (see #128).

### Changed
- Overhaul of the marker table structure and database build procedure (see #106, #111, #112).
- Use Byrska-Bishop 2022 instead of Auton 2015 for main frequency and Ae calculations (see #106, #112, #114, #130).
- Marker FASTA output now produces a single sequence per locus, even when multiple markers at a locus are defined (#124).


## [0.8.2] 2022-11-21

### Fixed
- Test suite execution bug.

## [0.8.1] 2022-11-21

### Fixed
- Bundling of files used for automated test suite.


## [0.8] 2022-11-21

### Added
- 25 marker definitions from Sun et al 2020 (see #87).
- 23 marker definitions from Jin et al 2020 (see #89).
- 20 marker definitions from Kureshi et al 2020 (see #90).
- GRCh38 offsets to `microhapdb marker --format=offsets` (see #92).
- 59 marker definitions from Wu et al 2021 (see #96).
- 22 marker definitions from Fan et al 2022 (see #97).
- New sphinx-based documentation suite hosted on ReadTheDocs.org (see #101).

### Changed
- Minor improvements to database build and corresponding documentation (see #88).
- Major overhaul of the Python API (see #99).

### Fixed
- Bug with dereferencing IDs and ID cross-references (see #99, #100).


## [0.7] 2022-01-20

### Added
- Markers and population frequency data (1000 Genomes and proprietary) for the mMHseq 90-plex panel (see #68).
- New `--extend-mode` option when displaying markers in `fasta` or `detail` mode (see #72).
- New options for preparing MicroHapDB data for import into MicroHapulator (see #80).
- New option for preparing MicroHapDB frequency data for import into probgen tools like LRMix Studio or EuroForMix (see #82).

### Changed
- Support for pandas>=1.2 added, support for Python 3.6 dropped (see #74).
- The `--marker` and `--population` flags for `microhapdb frequency` now support multiple arguments (see #81).

### Fixed
- Error message when detail view is requested for a marker without frequency information (see #67).
- Problematic haplotype frequencies for 4 mMHseq markers including SNPs whose rsIDs are not present in the original 1KGP data (see #73).


## [0.6] 2020-06-25

### Added
- Added F_ST (fixation index) to the markers table (see #59, #62, #63).
- Added the ability to swap out population-specific Ae values for the default 26-population average (see #63).
- Support for GRCh37 (see #64).

### Changed
- Various updates to documentation (see #60).
- Cosmetic updates to "marker --detail" (see #61).
- Refactoring of database build: cleanup using Pandas, reorganization for GRCh37 support (see #62, #64).


## [0.5] 2020-02-13

### Added
- Added 10 marker definitions from Voskoboinik et al 2018 (see #47).
- Added 118 marker definitions from de la Puente, Phillips, et al 2020 (see #49, #54).
- Added microhap frequency estimates for 26 global populations to all markers from the 1000 Genomes Project Phase 3 dat (see #55).
- Added I_n (informativeness for assignment) scores to each marker (see #56).

### Fixed
- Removed duplicate entry for marker `mh11CP-004` (see #46).
- Corrected marker name from `mh05KK-058` to `mh15KK-058` (see #46).
- Corrected erroneous rsID assignments for van der Gaag 2018 markers (see #46).

### Changed
- 26 global populations from the 1000 Genomes Project migrated from ALFRED data source to a dedicated 1KGP data source (see #45).
- Dropped Travis CI configuration (see #47).
- Changed default delta from 25 to 10 and default minlength from 250 to 80 (see #52).


## [0.4.3] 2019-11-05

### Fixed
- Squashed a bug with package data not being installed correctly.


## [0.4.2] 2019-11-05

### Changed
- Some aspects of panel design and definition migrated from MicroHapulator to MicroHapDB (see #34, #35).
- Minor changes to `TargetAmplicon` class to support MicroHapulator (see #36).


## [0.4.1] 2019-10-22

### Fixed
- Corrected a copy/paste error (see #33).


## [0.4] 2019-10-22

### Added
- Marker definitions and allele frequencies for a 40 microhap panel presented at ISFG 2019 (see #25).
- Marker definitions and allele frequencies for 27 microhaps published in a 2015 *Legal Medicine* paper (see #31).
- Marker definitions and allele frequencies for 11 microhaps published in a 2019 *FSI: Genetics* paper (see #32).

### Changed
- Replaced all references to `locus` and `loci` with `marker` and `markers` in the main code base (see #23).
- Huge overhaul to the database build procedure, the CLI, and the Python API (see #27).


## [0.3] 2019-05-02

### Added
- Marker definitions and allele frequency data for 16 microhaps from a 2018 *FSI: Genetics* paper (see #14).
- Function to compute SNP positions from microhap locus ID (see #18).
- A function to compute standard internal MicroHapDB ID for a single label or list of labels (see 9ec1e93735);
  any combination of ALFRED, LOVD, and MicroHapDB identifiers are valid input.

### Changed
- Replaced pip/PyPI installation instructions with bioconda installation instructions (see #14).

### Fixed
- Corrected DB build and coordinates (see #15).



## [0.2] 2018-12-06

### Changed
- New command-line interface.
- New Python API.
- New database based on an updated table of 198 microhaplotype loci across 96 populations from ALFRED.


## [0.1.2] 2018-10-10

Fixed a bug with testing the installed package.


## [0.1.1] 2018-10-10

Fixed a bug with distributing non-code files.


## [0.1] 2018-10-10

Initial release!
