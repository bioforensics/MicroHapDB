## help:      print this help message and exit
help: Makefile
	@echo ''
	@sed -n 's/^## //p' Makefile
	@echo ''

## test:      execute the automated test suite
test:
	pytest --cov=microhapdb --cov-report=term --cov-report=xml --doctest-modules microhapdb/cli/*.py microhapdb/marker.py microhapdb/population.py microhapdb/tests/test_*.py

## devdeps:   install development dependencies
devdeps:
	pip install --upgrade pip setuptools
	pip install wheel twine
	pip install black==22.10 'pytest>=5.0' pytest-cov

## clean:     remove development artifacts
clean:
	rm -rf __pycache__/ microhapdb/__pycache__/ microhapdb/*/__pycache__ build/ dist/ *.egg-info/ dbbuild/.snakemake

## style:     check code style
style:
	black --line-length=99 --check *.py microhapdb/*.py microhapdb/*/*.py

## format:    autoformat code with Black
format:
	black --line-length=99 *.py microhapdb/*.py microhapdb/*/*.py
