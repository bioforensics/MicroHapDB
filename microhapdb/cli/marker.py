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

from argparse import RawDescriptionHelpFormatter
import microhapdb
from microhapdb import Marker
import pandas as pd
import sys
from textwrap import dedent
from warnings import warn


def main(args):
    if args.ae_pop:
        microhapdb.set_ae_population(popid=args.ae_pop)
    if args.GRCh37:
        microhapdb.set_reference(37)
    markerids = resolve_panel(args.panel) if args.panel else args.id
    result = apply_filters(markerids, args.region, args.query)
    display(
        result,
        args.format,
        delta=args.delta,
        minlen=args.min_length,
        extend_mode=args.extend_mode,
        trunc=args.trunc,
        refr37=args.GRCh37,
    )
    if args.ae_pop:
        microhapdb.set_ae_population(popid=None)  # Reset
    if args.GRCh37:
        microhapdb.set_reference(38)  # Reset


def resolve_panel(panel):
    markerids = list()
    if hasattr(microhapdb.panel, args.panel):
        func = getattr(microhapdb.panel, args.panel)
        markerids.extend(func())
    else:
        with open(args.panel, "r") as fh:
            markerids = fh.read().strip().split()
    return markerids


def apply_filters(markerids=None, region=None, query=None):
    result = microhapdb.markers
    if region:
        result = Marker.table_from_region(region)
    if query:
        result = result.query(query, engine="python")
    if markerids:
        markerids = Marker.standardize_ids(markerids)
        result = result[result.Name.isin(markerids)]
    return result


def display(result, view_format, delta=10, minlen=80, extend_mode=0, trunc=True, refr37=False):
    if view_format == "table":
        if trunc:
            print(result.to_string(index=False))
        else:
            result.to_csv(sys.stdout, sep="\t", index=False)
    else:
        markers = list(
            Marker.objectify(result, delta=delta, minlen=minlen, extendmode=extend_mode)
        )
        if view_format == "detail":
            for marker in markers:
                print(marker.detail)
        elif view_format == "fasta":
            for marker in markers:
                print(marker.fasta)
        elif view_format == "offsets":
            refr = "Hg37" if refr37 else "Hg38"
            table = pd.concat([marker.definition for marker in markers])
            table = table.rename(columns={"ChromOffset": f"Offset{refr}"})
            table.to_csv(sys.stdout, sep="\t", index=False)
        else:
            raise ValueError(f'unsupported view format "{args.format}"')


def subparser(subparsers):
    desc = microhapdb.cli.bubbletext + "\nRetrieve marker records by identifier or query"
    epilog = """\
    Examples::

        microhapdb marker mh01NK-001
        microhapdb marker --format=fasta mh13KK-218 mh04CP-002 mh02AT-05
        microhapdb marker --format=fasta --panel mypanel.txt
        microhapdb marker --format=detail --min-length=125 --extend-mode=3 MHDBM-dc55cd9e
        microhapdb marker --region=chr18:1-25000000 --GRCh37
        microhapdb marker --query='Source == "ALFRED"' --ae-pop CEU
        microhapdb marker --query='Name.str.contains("PK")'
    """
    epilog = dedent(epilog)
    subparser = subparsers.add_parser(
        "marker",
        description=desc,
        epilog=epilog,
        formatter_class=RawDescriptionHelpFormatter,
    )
    subparser.add_argument(
        "--format", choices=["table", "detail", "fasta", "offsets"], default="table"
    )
    subparser.add_argument(
        "--ae-pop",
        metavar="POP",
        help="specify the 1000 Genomes population from which to report "
        'effective number of alleles in the "Ae" column; by default, the Ae value averaged over '
        "all 26 1KGP populations is reported",
    )
    subparser.add_argument(
        "--GRCh37",
        action="store_true",
        help="use coordinates from the GRCh37 reference "
        "assembly; by default, the GRCh38 reference is used",
    )
    subparser.add_argument(
        "--delta",
        metavar="D",
        type=int,
        default=10,
        help="extend D nucleotides beyond the "
        "marker extent when computing amplicon boundaries (detail and fasta format only); by "
        "default D=10",
    )
    subparser.add_argument(
        "--min-length",
        metavar="L",
        type=int,
        default=80,
        help="minimum amplicon length (detail " "and fasta format only); by default L=80",
    )
    subparser.add_argument(
        "--extend-mode",
        metavar="E",
        type=str_to_extend_mode,
        default="symmetric",
        help="specify how coordinates will be adjusted if extension is required to satisfy the "
        "minimum amplicon length; use `5` to extend the 5' end, `3` to extend the 3' end, or "
        "`symmetric` to extend both ends equally; by default, symmetric mode is used",
    )
    subparser.add_argument(
        "--notrunc",
        dest="trunc",
        action="store_false",
        default=True,
        help="disable truncation of tabular results",
    )
    subparser.add_argument(
        "-p",
        "--panel",
        metavar="FILE",
        help="file containing a list of marker names/identifiers," " one per line",
    )
    subparser.add_argument(
        "-r",
        "--region",
        metavar="RGN",
        help="restrict results to the " "specified genomic region; format chrX:YYYY-ZZZZZ",
    )
    subparser.add_argument(
        "--query", metavar="STRING", help="Retrieve records using a Pandas-style query"
    )
    subparser.add_argument("id", nargs="*", help="marker identifier")


def str_to_extend_mode(value):
    value = str(value)
    if value not in ("5", "3", "symmetric"):
        raise ValueError("extend mode must be `5`, `3`, or `symmetric`")
    if value == "5":
        return -1
    elif value == "3":
        return 1
    else:
        return 0
