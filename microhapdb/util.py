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
