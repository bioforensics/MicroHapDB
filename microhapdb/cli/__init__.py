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

from . import lookup, marker, population, frequency, summarize
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import microhapdb
from pyfaidx import Fasta as FastaIdx
from subprocess import run
import sys
from urllib.request import urlretrieve


subparser_funcs = {
    "lookup": lookup.subparser,
    "marker": marker.subparser,
    "population": population.subparser,
    "frequency": frequency.subparser,
    "summarize": summarize.subparser,
}

mains = {
    "lookup": lookup.main,
    "marker": marker.main,
    "population": population.main,
    "frequency": frequency.main,
    "summarize": summarize.main,
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
        raise SystemExit()
    if args.download:  # pragma: no cover
        download_hg38()
        raise SystemExit()
    assert args.cmd in mains
    mainmethod = mains[args.cmd]
    mainmethod(args)


def get_parser():
    cli = ArgumentParser(description=bubbletext, formatter_class=RawDescriptionHelpFormatter)
    cli._positionals.title = "Subcommands"
    cli._optionals.title = "Configuration"
    cli.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"MicroHapDB v{microhapdb.__version__}",
    )
    cli.add_argument(
        "-f", "--files", action="store_true", help="print data table filenames and exit"
    )
    cli.add_argument("--download", action="store_true", help="download the GRCh38 genome and exit")
    subcommandstr = ", ".join(sorted(subparser_funcs.keys()))
    subparsers = cli.add_subparsers(dest="cmd", metavar="cmd", help=subcommandstr)
    for func in subparser_funcs.values():
        func(subparsers)
    return cli


def print_files():
    tables = ("marker", "marker-aes", "population")
    for table in tables:
        print(microhapdb.data_file(f"{table}.csv"))
    print(microhapdb.data_file("frequency.csv.gz"))


def download_hg38():  # pragma: no cover
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    path = microhapdb.data_file("hg38.fasta")
    pathgz = f"{path}.gz"
    urlretrieve(url, pathgz)
    run(["gunzip", pathgz])
    hg38 = FastaIdx(path)
    markerseq = hg38["chr13"][53486574:53486837]  # Ensure index is built
