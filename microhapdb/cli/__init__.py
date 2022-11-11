# -*- coding: utf-8 -*-
#
# -------------------------------------------------------------------------------------------------
# Copyright (c) 2019, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from . import lookup, marker, population, frequency
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import microhapdb
from pkg_resources import resource_filename
import sys
import textwrap


subparser_funcs = {
    "lookup": lookup.subparser,
    "marker": marker.subparser,
    "population": population.subparser,
    "frequency": frequency.subparser,
}

mains = {
    "lookup": lookup.main,
    "marker": marker.main,
    "population": population.main,
    "frequency": frequency.main,
}

bubbletext = r"""
≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
 __  __ _            _  _           ___  ___
|  \/  (_)__ _ _ ___| || |__ _ _ __|   \| _ )
| |\/| | / _| '_/ _ \ __ / _` | '_ \ |) | _ \
|_|  |_|_\__|_| \___/_||_\__,_| .__/___/|___/
                              |_|
≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
"""


def main(args=None):
    if args is None:  # pragma: no cover
        args = get_parser().parse_args()
    try:  # If no arguments are provided, invoke --help mode
        if set([getattr(args, key) for key in vars(args)]) == set([False, None]):
            get_parser().parse_args(["-h"])
    except TypeError:
        pass
    if args.files:
        print_files()
        raise SystemExit
    assert args.cmd in mains
    mainmethod = mains[args.cmd]
    mainmethod(args)


def get_parser():
    cli = ArgumentParser(description=bubbletext, formatter_class=RawDescriptionHelpFormatter)
    cli._optionals.title = "Configuration"
    cli.add_argument(
        "-v",
        "--version",
        action="version",
        version="MicroHapDB v{}".format(microhapdb.__version__),
    )
    cli.add_argument(
        "-f", "--files", action="store_true", help="print data table filenames and exit"
    )
    subcommandstr = ", ".join(sorted(subparser_funcs.keys()))
    subparsers = cli.add_subparsers(dest="cmd", metavar="cmd", help=subcommandstr)
    for func in subparser_funcs.values():
        func(subparsers)
    return cli


def print_files():
    tables = ("marker", "population", "frequency", "variantmap", "idmap", "sequence", "indels")
    for table in tables:
        print(resource_filename("microhapdb", f"{table}.tsv"))
