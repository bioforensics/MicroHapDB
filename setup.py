# -------------------------------------------------------------------------------------------------
# Copyright (c) 2018, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from setuptools import setup
import versioneer


desc = "Portable database of microhaplotype marker and allele frequency data"
with open("README.md", "r") as infile:
    longdesc = infile.read()

setup(
    name="microhapdb",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description=desc,
    long_description=longdesc,
    long_description_content_type="text/markdown",
    url="https://github.com/bioforensics/MicroHapDB",
    author="Daniel Standage",
    author_email="daniel.standage@nbacc.dhs.gov",
    packages=["microhapdb", "microhapdb.cli", "microhapdb.tests"],
    package_data={"microhapdb": ["microhapdb/data/*", "microhapdb/data/tests/*"]},
    include_package_data=True,
    install_requires=["pandas>=1.2", "pyfaidx>=0.7"],
    entry_points={"console_scripts": ["microhapdb = microhapdb.cli:main"]},
    classifiers=[
        "Environment :: Console",
        "Framework :: IPython",
        "Framework :: Jupyter",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Legal Industry",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    zip_safe=True,
)
