#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse
import microhapdb
from microhapdb.retrieve import fetch_by_id, fetch_by_query, fetch_by_region
import sys
import textwrap


def query_mode(table, querystr):
    """Retrieve data with a Pandas query.

    Specify the table to query (variant, locus, population, allele) and provide
    a Pandas-style query.

    >>> for result in query_mode('locus', 'Chrom == "chr8"'): print(result)
                  ID Reference Chrom     Start       End  Source
    77   MHDBL000078    GRCh38  chr8  11738320  11738460  ALFRED
    140  MHDBL000141    GRCh38  chr8   3659270   3659482  ALFRED
    >>> querystr = 'Locus == "MHDBL000177" and Population == "MHDBP000092"'
    >>> for result in query_mode('allele', querystr): print(result)
                 Locus   Population   Allele  Frequency
    71439  MHDBL000177  MHDBP000092  T,C,A,A      0.288
    71440  MHDBL000177  MHDBP000092  T,C,G,A      0.025
    71441  MHDBL000177  MHDBP000092  T,C,G,G      0.071
    71442  MHDBL000177  MHDBP000092  T,G,A,A      0.268
    71443  MHDBL000177  MHDBP000092  T,G,A,G      0.005
    71444  MHDBL000177  MHDBP000092  T,G,G,A      0.025
    71445  MHDBL000177  MHDBP000092  T,G,G,G      0.318
    71446  MHDBL000177  MHDBP000092  G,C,A,A      0.000
    71447  MHDBL000177  MHDBP000092  G,C,A,G      0.000
    71448  MHDBL000177  MHDBP000092  G,C,G,A      0.000
    71449  MHDBL000177  MHDBP000092  G,C,G,G      0.000
    71450  MHDBL000177  MHDBP000092  G,G,A,A      0.000
    71451  MHDBL000177  MHDBP000092  G,G,G,A      0.000
    71452  MHDBL000177  MHDBP000092  G,G,G,G      0.000
    """
    for result in fetch_by_query(table, querystr):
        yield result


def id_mode(idstr):
    """Retrieve data using internal or external IDs, names, or labels.

    >>> for result in id_mode('SI605775E'): print(result)
                      ID Reference  Chrom  Position Alleles    Source
    9961  MHDBV000009962    GRCh38  chr13  50313423     C,T  dbSNP151
    >>> for result in id_mode('rs690302'): print(result)
                       ID Reference  Chrom  Position  Alleles    Source
    17205  MHDBV000017206    GRCh38  chr18   8892896  A,C,G,T  dbSNP151
    >>> for result in id_mode('mh19CP-007'): print(result)
                 ID Reference  Chrom     Start       End  Source
    66  MHDBL000067    GRCh38  chr19  14310740  14310781  ALFRED
    >>> for result in id_mode('SA004109O'): print(result)
                 ID       Name  Source
    76  MHDBP000077  Colombian  ALFRED
    """
    for result in fetch_by_id(idstr):
        yield result


def region_mode(region, table=None):
    """Retrieve microhap loci and proximal variants with range queries.

    Use queries of the format "chrX" or "chrX:YYYY-ZZZZ" to retrieve loci or
    variants (or both) from the specified genomic region.

    MicroHapDB includes not only the dbSNP variants that define each
    microhaplotype locus, but also the other variants within its extent and the
    flanking nucleotides.

    >>> for result in region_mode('chr13', table='locus'): print(result)
                  ID Reference  Chrom      Start        End  Source
    44   MHDBL000045    GRCh38  chr13   24343963   24343994  ALFRED
    50   MHDBL000051    GRCh38  chr13   23191402   23191542  ALFRED
    51   MHDBL000052    GRCh38  chr13   46291795   46291987  ALFRED
    52   MHDBL000053    GRCh38  chr13   53486692   53486837  ALFRED
    68   MHDBL000069    GRCh38  chr13  110154352  110154505  ALFRED
    69   MHDBL000070    GRCh38  chr13   66138600   66138696  ALFRED
    70   MHDBL000071    GRCh38  chr13   94894396   94894513  ALFRED
    116  MHDBL000117    GRCh38  chr13   50313424   50313589  ALFRED
    >>> for result in region_mode('chr18:24557400-24557450'): print(result)
                 ID Reference  Chrom     Start       End  Source
    34  MHDBL000035    GRCh38  chr18  24557355  24557490  ALFRED
                       ID Reference  Chrom  Position Alleles    Source
    17349  MHDBV000017350    GRCh38  chr18  24557400     C,T  dbSNP151
    17350  MHDBV000017351    GRCh38  chr18  24557402     A,C  dbSNP151
    17351  MHDBV000017352    GRCh38  chr18  24557414     A,T  dbSNP151
    17352  MHDBV000017353    GRCh38  chr18  24557416     A,G  dbSNP151
    17353  MHDBV000017354    GRCh38  chr18  24557431     A,G  dbSNP151
    17354  MHDBV000017355    GRCh38  chr18  24557443     A,G  dbSNP151
    17355  MHDBV000017356    GRCh38  chr18  24557447   C,G,T  dbSNP151
    17356  MHDBV000017357    GRCh38  chr18  24557448     A,G  dbSNP151
    """
    for result in fetch_by_region(region, table):
        yield result


def print_files():
    """Print the location of the installed data files on the system.

    This is helpful if the user wants to data science the TSV tables directly.
    """
    tables = (
        'allele', 'locus', 'population', 'variant', 'variantmap', 'idmap'
    )
    for table in tables:
        print(microhapdb.data_file(table + '.tsv'))


def get_parser():
    """Construct an argument parser for the command-line interface."""
    bubbletext = r'''
≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
 __  __ _            _  _           ___  ___
|  \/  (_)__ _ _ ___| || |__ _ _ __|   \| _ )
| |\/| | / _| '_/ _ \ __ / _` | '_ \ |) | _ \
|_|  |_|_\__|_| \___/_||_\__,_| .__/___/|___/
                              |_|
≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠≠
'''
    epilog = """\
    Examples::

        microhapdb --id mh12KK-043
        microhapdb --id rs10815466
        microhapdb --id SA000936S
        microhapdb --table population
        microhapdb --table locus --region chr7
        microhapdb --table variant --region chr19:4000000-5000000
        microhapdb --table allele --query 'Locus == "MHDBL000047" and Allele == "C,T,G" and Population == "MHDBP000006"'
    """
    epilog = textwrap.dedent(epilog)
    cli = argparse.ArgumentParser(
        description=bubbletext, epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    cli._optionals.title = 'Configuration'
    cli.add_argument(
        '-v', '--version', action='version',
        version='MicroHapDB v{}'.format(microhapdb.__version__)
    )
    cli.add_argument(
        '-f', '--files', action='store_true', help='print data table '
        'filenames and exit'
    )
    cli.add_argument(
        '-t', '--table', choices=['variant', 'locus', 'population', 'allele'],
        metavar='TBL', help='restrict results to the specified data table; '
        'must be one of "variant", "locus", "population", or "allele"'
    )
    cli.add_argument(
        '-r', '--region', metavar='RGN', help='restrict results to the '
        'specified genomic region; format chrX:YYYY-ZZZZZ'
    )
    cli.add_argument(
        '--id', metavar='ID', help='query data tables using a dbSNP ID, an '
        'ALFRED ID/name, or an internal MicroHapDB ID'
    )
    cli.add_argument(
        '-q', '--query', metavar='QRY', help='Invoke a Pandas-style query; '
        'must specify table to query with the `-t|--table` flag'
    )

    return cli


def main(args=None):
    """MicroHapDB main method."""
    if args is None:  # pragma: no cover
        args = get_parser().parse_args()

    # If no arguments are provided, invoke --help mode.
    if set([getattr(args, key) for key in vars(args)]) == set([False, None]):
        get_parser().parse_args(['-h'])

    if args.files:
        print_files()
        raise SystemExit(0)

    # Who's a good doggy?
    retriever = None
    if args.query:
        for attr in ('region', 'id'):
            if getattr(args, attr):
                msg = 'ignoring "{}" parameter in "query" mode'.format(attr)
                print('[MicroHapDB] WARNING:', msg, file=sys.stderr)
        retriever = query_mode(args.table, args.query)
    elif args.id:
        if args.table:
            message = 'ignoring "table" parameter in "id" mode'
            print('[MicroHapDB] WARNING:', message, file=sys.stderr)
        retriever = id_mode(args.id)
    elif args.region:
        retriever = region_mode(args.region, args.table)
    else:
        retriever = [microhapdb.tables[args.table]]

    for result in retriever:
        print(result.to_string())
