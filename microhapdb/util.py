#!/usr/bin/env python3
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import os
from pkg_resources import resource_filename


def data_file(path):
    """Return full path to a MicroHapDB data file

    Data files are installed *in situ* along with MicroHapDB's Python code.
    This helper function is used internally to determine the full path from the
    package's installation point.
    """
    pathparts = path.split('/')
    relpath = os.path.join('data', *pathparts)
    return resource_filename('microhapdb', relpath)


def parse_regionstr(regionstr):
    '''Retrieve chromosome name and coordinates from a region string

    Region string is expected to be in one of the two following formats: 'chr3'
    or 'chr3:1000000-5000000'.

    >>> parse_regionstr('chr12')
    ('chr12', None, None)
    >>> parse_regionstr('chr12:345-678')
    ('chr12', 345, 678)
    '''
    chrom, start, end = None, None, None
    if ':' in regionstr:
        chrom, rng = regionstr.split(':')
        if rng.count('-') != 1:
            raise ValueError('cannot parse region "{}"'.format(regionstr))
        startstr, endstr = rng.split('-')
        start, end = int(startstr), int(endstr)
    else:
        chrom = regionstr
    return chrom, start, end
