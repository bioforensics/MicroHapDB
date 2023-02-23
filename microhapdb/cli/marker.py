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
    markerids = resolve_panel(args.panel) if args.panel else args.id
    result = apply_filters(markerids, args.region, args.query)
    if len(result) > 0:
        display(
            result,
            args.format,
            columns=args.columns,
            delta=args.delta,
            minlen=args.min_length,
            extend_mode=args.extend_mode,
            trunc=args.trunc,
        )
    if args.ae_pop:
        # Reset
        microhapdb.set_ae_population(popid="1KGP")


def resolve_panel(panel):
    markerids = list()
    if hasattr(microhapdb.panel, panel):
        func = getattr(microhapdb.panel, panel)
        markerids.extend(func())
    else:
        with open(panel, "r") as fh:
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


def display(
    result,
    view_format,
    columns="nxcse",
    delta=10,
    minlen=80,
    extend_mode=0,
    trunc=True,
):
    if view_format == "table":
        result = subset_result(result, columns)
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
            table = pd.concat([marker.definition for marker in markers])
            table = table.rename(columns={"ChromOffset": f"OffsetHg38"})
            table.to_csv(sys.stdout, sep="\t", index=False)
        else:
            raise ValueError(f'unsupported view format "{view_format}"')


def subset_result(result, columns):
    fmt = {
        "n": "NumVars",
        "x": "Extent",
        "c": "Chrom",
        "s": "Start",
        "e": "End",
        "p": "Positions",
        "q": "Positions37",
        "r": "RSIDs",
        "a": "Ae",
    }
    for code in columns:
        if code not in fmt:
            raise ValueError(f"unsupported format code '{code}'")
    cols = ["Name"] + [fmt[code] for code in columns] + ["Source"]
    return result[cols]


def subparser(subparsers):
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
        description="Retrieve marker records by identifier or query",
        epilog=epilog,
        formatter_class=RawDescriptionHelpFormatter,
    )
    retrieval = subparser.add_argument_group(
        "Data Retrieval", "Configure how marker records are retrieved from the database."
    )
    retrieval.add_argument(
        "--ae-pop",
        metavar="POP",
        type=str,
        help="specify the 1000 Genomes population from which to report "
        'effective number of alleles in the "Ae" column; by default, the Ae value averaged over '
        "all 26 1KGP populations is reported",
    )
    retrieval.add_argument(
        "--panel",
        metavar="FILE",
        help="file containing a list of marker names/identifiers," " one per line",
    )
    retrieval.add_argument(
        "--region",
        metavar="RGN",
        help="restrict results to the " "specified genomic region; format chrX:YYYY-ZZZZZ",
    )
    retrieval.add_argument(
        "--query", metavar="QRY", help="Retrieve records using a Pandas-style query"
    )
    formatting = subparser.add_argument_group(
        "Formatting",
        "Configure how results are formatted. Some formats include information for a 'target sequence' for each marker, representing what would be targeted by e.g. hybridization capture probes or PCR primers for amplicon sequencing. MicroHapDB computes the endpoints of these target sequences by extending `--delta=D` nucleotides beyond the first and last SNPs defining the marker, and then—if needed—extending further until `--min-length=L` is satisfied. Configuration of these and related parameters is described below.",
    )
    formatting.add_argument(
        "--format", choices=["table", "detail", "fasta", "offsets"], default="table"
    )
    formatting.add_argument(
        "--columns",
        metavar="C",
        default="nxcsea",
        help="string of column codes indicating which fields to include in tabular output; n=NumVars x=Extent c=Chrom s=Start e=End p=Positions q=Positions37 r=RSIDs a=Ae; by default C=nxcsea",
    )
    formatting.add_argument(
        "--delta",
        metavar="D",
        type=int,
        default=10,
        help="extend D nucleotides beyond the marker extent when computing target sequence boundaries; by default D=10",
    )
    formatting.add_argument(
        "--min-length",
        metavar="L",
        type=int,
        default=80,
        help="minimum length of the target sequence; by default L=80",
    )
    formatting.add_argument(
        "--extend-mode",
        metavar="E",
        type=str_to_extend_mode,
        default="symmetric",
        help="specify how the target sequence will be extended to satisfy the minimum length criterion; use `5` to extend only the 5' end, `3` to extend only the 3' end, or `symmetric` to extend both ends equally; by default, symmetric mode is used",
    )
    formatting.add_argument(
        "--notrunc",
        dest="trunc",
        action="store_false",
        default=True,
        help="disable truncation of tabular results",
    )
    subparser.add_argument("id", nargs="*", help="one or more marker identifiers")
    subparser._positionals.title = "Required Arguments"
    subparser._optionals.title = "Options"


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
