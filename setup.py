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
    'Resource for working with microhaploptype data from the ALFRED database '
    '(https://alfred.med.yale.edu).'
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
    url='https://github.com/bioforensics/microhapdb',
    author='Daniel Standage',
    author_email='daniel.standage@nbacc.dhs.gov',
    packages=['microhapdb'],
    package_data={
        'microhapdb': ['microhapdb/data/*', 'microhapdb/data/raw/*']
    },
    include_package_data=True,
    install_requires=['pandas', 'pytest'],
    entry_points={
        'console_scripts': ['microhapdb = microhapdb.cli:main']
    },
    classifiers=[
        'Environment :: Console',
        'Framework :: IPython',
        'Framework :: Jupyter',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    zip_safe=True,
)
