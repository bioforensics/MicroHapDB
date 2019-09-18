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

    Specify the table to query (variant, marker, population, allele) and
    provide a Pandas-style query.

    >>> for result in query_mode('marker', 'Chrom == "chr8"'): print(result)
                  ID Reference Chrom     Start       End   AvgAe  Source
    202  MHDBM000203    GRCh38  chr8   1194352   1194372  2.6273    LOVD
    203  MHDBM000204    GRCh38  chr8   3659269   3659482  3.9708  ALFRED
    204  MHDBM000205    GRCh38  chr8  11738319  11738460  2.2865  ALFRED
    >>> querystr = 'Marker == "MHDBM000177" and Population == "MHDBP000092"'
    >>> for result in query_mode('allele', querystr): print(result)
                Marker   Population Allele  Frequency
    66015  MHDBM000177  MHDBP000092    C,C      0.661
    66016  MHDBM000177  MHDBP000092    C,T      0.339
    66017  MHDBM000177  MHDBP000092    T,C      0.000
    66018  MHDBM000177  MHDBP000092    T,T      0.000
    """
    for result in fetch_by_query(table, querystr):
        yield result


def id_mode(idstr):
    """Retrieve data using internal or external IDs, names, or labels.

    >>> for result in id_mode('SI605775E'): print(result)
                       ID Reference  Chrom  Position Alleles    Source
    10491  MHDBV000010492    GRCh38  chr13  50313423     C,T  dbSNP151
    >>> for result in id_mode('rs690302'): print(result)
                       ID Reference  Chrom  Position  Alleles    Source
    19313  MHDBV000019314    GRCh38  chr18   8892896  A,C,G,T  dbSNP151
    >>> for result in id_mode('mh19CP-007'): print(result)
                  ID Reference  Chrom     Start       End   AvgAe  Source
    109  MHDBM000110    GRCh38  chr19  14310739  14310781  3.0813  ALFRED
    >>> for result in id_mode('SA004109O'): print(result)
                 ID       Name  Source
    15  MHDBP000016  Colombian  ALFRED
    """
    for result in fetch_by_id(idstr):
        yield result


def region_mode(region, table=None):
    """Retrieve microhap markers and proximal variants with range queries.

    Use queries of the format "chrX" or "chrX:YYYY-ZZZZ" to retrieve markers or
    variants (or both) from the specified genomic region.

    MicroHapDB includes not only the dbSNP variants that define each
    microhaplotype marker, but also the other variants within its extent and
    the flanking nucleotides.

    >>> for result in region_mode('chr13', table='marker'): print(result)
                 ID Reference  Chrom      Start        End   AvgAe  Source
    55  MHDBM000056    GRCh38  chr13   23191401   23191542  3.6440  ALFRED
    56  MHDBM000057    GRCh38  chr13   24343962   24343994  3.0655  ALFRED
    57  MHDBM000058    GRCh38  chr13   46291794   46291987  4.0035  ALFRED
    58  MHDBM000059    GRCh38  chr13   50313423   50313589  2.3932  ALFRED
    59  MHDBM000060    GRCh38  chr13   53486691   53486837  6.0444  ALFRED
    60  MHDBM000061    GRCh38  chr13   66138599   66138696  3.4361  ALFRED
    61  MHDBM000062    GRCh38  chr13   94894395   94894513  2.0150  ALFRED
    62  MHDBM000063    GRCh38  chr13  110154351  110154505  3.7978  ALFRED
    >>> for result in region_mode('chr18:24557400-24557450'): print(result)
                  ID Reference  Chrom     Start       End  AvgAe  Source
    105  MHDBM000106    GRCh38  chr18  24557354  24557490  2.658  ALFRED
                       ID Reference  Chrom  Position Alleles    Source
    19457  MHDBV000019458    GRCh38  chr18  24557400     C,T  dbSNP151
    19458  MHDBV000019459    GRCh38  chr18  24557402     A,C  dbSNP151
    19459  MHDBV000019460    GRCh38  chr18  24557414     A,T  dbSNP151
    19460  MHDBV000019461    GRCh38  chr18  24557416     A,G  dbSNP151
    19461  MHDBV000019462    GRCh38  chr18  24557431     A,G  dbSNP151
    19462  MHDBV000019463    GRCh38  chr18  24557443     A,G  dbSNP151
    19463  MHDBV000019464    GRCh38  chr18  24557447   C,G,T  dbSNP151
    19464  MHDBV000019465    GRCh38  chr18  24557448     A,G  dbSNP151
    """
    for result in fetch_by_region(region, table):
        yield result


def print_files():
    """Print the location of the installed data files on the system.

    This is helpful if the user wants to data science the TSV tables directly.
    """
    tables = (
        'allele', 'marker', 'population', 'variant', 'variantmap', 'idmap'
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
        microhapdb --table marker --region chr7
        microhapdb --table variant --region chr19:4000000-5000000
        microhapdb --table allele --query 'Marker == "MHDBM000047" and Allele == "C,T,G" and Population == "MHDBP000006"'
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
        '-t', '--table', choices=['variant', 'marker', 'population', 'allele'],
        metavar='TBL', help='restrict results to the specified data table; '
        'must be one of "variant", "marker", "population", or "allele"'
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
