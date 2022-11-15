# Installation

## Overview

For best results, installation from [bioconda](https://bioconda.github.io/) is recommended.

```
conda create --name microhapdb -y python=3.9 microhapdb
conda activate microhapdb
microhapdb --help
```

To make sure the package installed correctly:

```
conda install -y pytest
pytest --pyargs microhapdb --doctest-modules
```

Installation with conda ensures that compatible versions of Python and the [Pandas](https://pandas.pydata.org) library are correctly installed.


## Development quick start

If you're setting up an environment to contribute to development of MicroHapDB's CLI or API, you'll want to skip the procedure outlined above and use the following instead.

```
conda create --name microhapdb -y python=3.9 
conda activate microhapdb
git clone https://github.com/bioforensics/MicroHapDB.git
cd MicroHapDB/
pip install -e .  # Install the package and its Python dependencies
make devdeps      # Install development packages
make devhooks     # Register pre-commit hooks for development
```


## Extending or rebuilding the database

Additional software dependencies and data files are required to add data to MicroHapDB or to reconstruct its contents from scratch.
See [https://github.com/bioforensics/MicroHapDB/tree/master/dbbuild](https://github.com/bioforensics/MicroHapDB/tree/master/dbbuild) for additional details.
