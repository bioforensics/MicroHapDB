# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## Unreleased

### Changed
- Replaced all references to `locus` and `loci` with `marker` and `markers` in the main code base (see #23).

### Fixed
- A bug in the dbbuild code that resulted in a variant off-by-one error (see #21).
- A bug in the `allele_positions()` code that didn't properly handle cases where known duplicate variants exist (see #21).


## [0.3] 2019-05-02

### Added
- Microhaps and limited allele frequency data from LOVD (see #14).
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
