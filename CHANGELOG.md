# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## Unreleased

### Added
- Added 10 marker definitions from Voskoboinik et al 2018 (see #47).
- Added 118 marker definitions from de la Puente, Phillips, et al 2020 (see #49, #54).

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
