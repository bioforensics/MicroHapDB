## help:      print this help message and exit
help: Makefile
	@echo ''
	@sed -n 's/^## //p' Makefile
	@echo ''

## test:      execute the automated test suite
test:
	pytest --cov=microhapdb --doctest-modules microhapdb/cli/*.py microhapdb/retrieve.py microhapdb/test_retrieve.py microhapdb/test_tables.py

## devdeps:   install development dependencies
devdeps:
	pip install --upgrade pip setuptools
	pip install wheel twine
	pip install pycodestyle 'pytest>=5.0' pytest-cov pytest-sugar

## clean:     remove development artifacts
clean:
	rm -rf __pycache__/ microhapdb/__pycache__/ microhapdb/*/__pycache__ build/ dist/ *.egg-info/ dbbuild/.snakemake

## style:     check code style against PEP8
style:
	pycodestyle --ignore=E501 microhapdb/*.py
