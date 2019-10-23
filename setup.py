#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from setuptools import setup
import versioneer


desc = (
    'Portable database of microhaplotype marker and allele frequency data'
)
with open('README.md', 'r') as infile:
    longdesc = infile.read()

setup(
    name='microhapdb',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description=desc,
    long_description=longdesc,
    long_description_content_type='text/markdown',
    url='https://github.com/bioforensics/MicroHapDB',
    author='Daniel Standage',
    author_email='daniel.standage@nbacc.dhs.gov',
    packages=['microhapdb', 'microhapdb.cli', 'microhapdb.tests'],
    package_data={
        'microhapdb': ['microhapdb/data/*', 'microhapdb/data/tests/*']
    },
    include_package_data=True,
    install_requires=['pandas', 'pytest>=5.0'],
    entry_points={
        'console_scripts': ['microhapdb = microhapdb.cli:main']
    },
    classifiers=[
        'Environment :: Console',
        'Framework :: IPython',
        'Framework :: Jupyter',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Legal Industry',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    zip_safe=True,
)
