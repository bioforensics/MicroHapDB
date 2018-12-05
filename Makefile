## help:      print this help message and exit
help: Makefile
	@echo ''
	@sed -n 's/^## //p' Makefile
	@echo ''

## test:      execute the automated test suite
test:
	py.test --cov=microhapdb --doctest-modules microhapdb/test_*.py microhapdb/retrieve.py

## devdeps:   install development dependencies
devdeps:
	pip3 install --upgrade pip setuptools
	pip3 install wheel twine
	pip3 install pycodestyle pytest-cov pytest-sugar

## clean:     remove development artifacts
clean:
	rm -rf __pycache__/ microhapdb/__pycache__/ microhapdb/*/__pycache__ build/ dist/ *.egg-info/ dbbuild/.snakemake

## style:     check code style against PEP8
style:
	pycodestyle --ignore=E501 microhapdb/*.py dbbuild/*.py dbbuild/Snakefile
