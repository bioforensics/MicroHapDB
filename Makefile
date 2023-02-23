## help:      print this help message and exit
help: Makefile
	@echo ''
	@sed -n 's/^## //p' Makefile
	@echo ''

## test:      execute the automated test suite
test:
	pytest --cov=microhapdb --cov-report=term --cov-report=xml --doctest-modules --pyargs microhapdb

## devdeps:   install development dependencies
devdeps:
	pip install --upgrade pip setuptools
	pip install wheel twine
	pip install black==22.10 'pytest>=5.0' pytest-cov myst-parser sphinx sphinx-argparse

## clean:     remove development artifacts
clean:
	rm -rf __pycache__/ microhapdb/__pycache__/ microhapdb/*/__pycache__ build/ dist/ *.egg-info/ dbbuild/.snakemake

## style:     check code style
style:
	black --line-length=99 --check *.py microhapdb/*.py microhapdb/*/*.py

## format:    autoformat code with Black
format:
	black --line-length=99 *.py microhapdb/*.py microhapdb/*/*.py

## doc:       build HTML documentation
doc:
	sphinx-build -b html docs docs/_build/

## hg38:      download and index the GRCh38 reference genome
hg38:
	curl -L -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > microhapdb/data/hg38.fasta.gz
	gunzip microhapdb/data/hg38.fasta.gz
	faidx microhapdb/data/hg38.fasta chr13:53486575-53486837

## devhooks:  install development hooks
devhooks:
	echo 'set -eo pipefail' > .git/hooks/pre-commit
	echo 'make style' >> .git/hooks/pre-commit
	echo 'make doc' >> .git/hooks/pre-commit
	chmod 755 .git/hooks/pre-commit
