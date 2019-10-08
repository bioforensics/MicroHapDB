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
                  ID Reference Chrom     Start       End   AvgAe                            Source
    238  MHDBM000239    GRCh38  chr8   1194352   1194372  2.6273  doi:10.1016/j.fsigen.2018.05.008
    239  MHDBM000240    GRCh38  chr8   3659269   3659482  3.9708                            ALFRED
    240  MHDBM000241    GRCh38  chr8  11738319  11738460  2.2865                            ALFRED
    >>> querystr = 'Marker == "MHDBM000213" and Population == "MHDBP000093"'
    >>> for result in query_mode('allele', querystr): print(result)
                Marker   Population Allele  Frequency
    66015  MHDBM000213  MHDBP000093    C,C      0.661
    66016  MHDBM000213  MHDBP000093    C,T      0.339
    66017  MHDBM000213  MHDBP000093    T,C      0.000
    66018  MHDBM000213  MHDBP000093    T,T      0.000
    """
    for result in fetch_by_query(table, querystr):
        yield result


def id_mode(idstr):
    """Retrieve data using internal or external IDs, names, or labels.

    >>> for result in id_mode('SI605775E'): print(result)
                       ID Reference  Chrom  Position Alleles    Source
    10265  MHDBV000010266    GRCh38  chr13  50313423     C,T  dbSNP151
    >>> for result in id_mode('rs690302'): print(result)
                       ID Reference  Chrom  Position  Alleles    Source
    18914  MHDBV000018915    GRCh38  chr18   8892896  A,C,G,T  dbSNP151
    >>> for result in id_mode('mh19CP-007'): print(result)
                  ID Reference  Chrom     Start       End   AvgAe  Source
    130  MHDBM000131    GRCh38  chr19  14310739  14310781  3.0813  ALFRED
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
                 ID Reference  Chrom      Start        End   AvgAe         Source
    65  MHDBM000066    GRCh38  chr13   23191401   23191542  3.6440         ALFRED
    66  MHDBM000067    GRCh38  chr13   23191461   23191542  2.5124  ISFG2019:P597
    67  MHDBM000068    GRCh38  chr13   24343962   24343994  3.0655         ALFRED
    68  MHDBM000069    GRCh38  chr13   46291794   46291987  4.0035         ALFRED
    69  MHDBM000070    GRCh38  chr13   46291948   46291987  2.6806  ISFG2019:P597
    70  MHDBM000071    GRCh38  chr13   50313423   50313589  2.3932         ALFRED
    71  MHDBM000072    GRCh38  chr13   53486691   53486837  6.0444         ALFRED
    72  MHDBM000073    GRCh38  chr13   66138599   66138696  3.4361         ALFRED
    73  MHDBM000074    GRCh38  chr13   94894395   94894513  2.0150         ALFRED
    74  MHDBM000075    GRCh38  chr13  110154351  110154412  3.2411  ISFG2019:P597
    75  MHDBM000076    GRCh38  chr13  110154351  110154505  3.7978         ALFRED
    >>> for result in region_mode('chr18:24557400-24557450'): print(result)
                  ID Reference  Chrom     Start       End   AvgAe         Source
    124  MHDBM000125    GRCh38  chr18  24557354  24557490  2.6580         ALFRED
    125  MHDBM000126    GRCh38  chr18  24557431  24557490  3.1874  ISFG2019:P597
                       ID Reference  Chrom  Position Alleles    Source
    19058  MHDBV000019059    GRCh38  chr18  24557400     C,T  dbSNP151
    19059  MHDBV000019060    GRCh38  chr18  24557402     A,C  dbSNP151
    19060  MHDBV000019061    GRCh38  chr18  24557414     A,T  dbSNP151
    19061  MHDBV000019062    GRCh38  chr18  24557416     A,G  dbSNP151
    19062  MHDBV000019063    GRCh38  chr18  24557431     A,G  dbSNP151
    19063  MHDBV000019064    GRCh38  chr18  24557443     A,G  dbSNP151
    19064  MHDBV000019065    GRCh38  chr18  24557447   C,G,T  dbSNP151
    19065  MHDBV000019066    GRCh38  chr18  24557448     A,G  dbSNP151
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
