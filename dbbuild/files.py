# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import builtins
from gzip import open as gzopen
import os
import re
from shutil import copyfileobj
from sys import stdin, stdout
from urllib.request import urlretrieve


def sopen(filename, mode):
    """Smart file handler

    Determines whether to create a compressed or un-compressed file handle
    based on filename extension.
    """
    if mode not in ('r', 'w'):
        raise ValueError('invalid mode "{}"'.format(mode))
    if filename in ['-', None]:
        filehandle = stdin if mode == 'r' else stdout
        return filehandle
    openfunc = builtins.open
    if filename.endswith('.gz'):
        openfunc = gzopen
        mode += 't'
    return openfunc(filename, mode)


def download_and_compress(url, destfile):
    """Download an uncompressed data file and compress it locally."""
    assert destfile.endswith('.gz')
    tempfile = os.path.splitext(destfile)[0]
    urlretrieve(url, tempfile)
    with open(tempfile, 'rb') as infile, gzopen(destfile, 'wb') as outfile:
        copyfileobj(infile, outfile)
    os.unlink(tempfile)


def cat(outstream, infiles):
    """Python drop-in for the UNIX `cat` command"""
    for infile in infiles:
        with sopen(infile, 'r') as instream:
            copyfileobj(instream, outstream)


def dlfile(path):
    """Determine full relative file path for a downloaded file.

    Downloaded files are placed in the `downloads` folder. Directory separators
    are handled in an OS-independent manner.
    """
    pathparts = re.split(r'[/\\]', path)
    return os.path.join('downloads', *pathparts)


def tmpfile(path):
    """Determine full relative file path for a processed intermediate file.

    Intermediate files created during the build process are placed in the
    `downloads` folder. Directory separators are handled in an OS-independent
    manner.
    """
    pathparts = re.split(r'[/\\]', path)
    return os.path.join('intermediate', *pathparts)
