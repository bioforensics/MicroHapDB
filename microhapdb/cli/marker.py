# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from argparse import RawDescriptionHelpFormatter
import microhapdb
from microhapdb.marker import print_detail, print_fasta, print_offsets
from textwrap import dedent
import sys


def str_to_extend_mode(value):
    value = str(value)
    if value not in ('5', '3', 'symmetric'):
        raise ValueError('extend mode must be `5`, `3`, or `symmetric`')
    if value == '5':
        return -1
    elif value == '3':
        return 1
    else:
        return 0


def subparser(subparsers):
    desc = microhapdb.cli.bubbletext + '\nRetrieve marker records by identifier or query'
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
        'marker', description=desc, epilog=epilog, formatter_class=RawDescriptionHelpFormatter,
    )
    subparser.add_argument('--format', choices=['table', 'detail', 'fasta', 'offsets'], default='table')
    subparser.add_argument(
        '--ae-pop', metavar='POP', help='specify the 1000 Genomes population from which to report '
        'effective number of alleles in the "Ae" column; by default, the Ae value averaged over '
        'all 26 1KGP populations is reported'
    )
    subparser.add_argument(
        '--GRCh37', action='store_true', help='use coordinates from the GRCh37 reference '
        'assembly; by default, the GRCh38 reference is used'
    )
    subparser.add_argument(
        '--delta', metavar='D', type=int, default=10, help='extend D nucleotides beyond the '
        'marker extent when computing amplicon boundaries (detail and fasta format only); by '
        'default D=10'
    )
    subparser.add_argument(
        '--min-length', metavar='L', type=int, default=80, help='minimum amplicon length (detail '
        'and fasta format only); by default L=80'
    )
    subparser.add_argument(
        '--extend-mode', metavar='E', type=str_to_extend_mode, default='symmetric',
        help="specify how coordinates will be adjusted if extension is required to satisfy the "
        "minimum amplicon length; use `5` to extend the 5' end, `3` to extend the 3' end, or "
        "`symmetric` to extend both ends equally; by default, symmetric mode is used"
    )
    subparser.add_argument(
        '--notrunc', dest='trunc', action='store_false', default=True,
        help='disable truncation of tabular results'
    )
    subparser.add_argument(
        '-p', '--panel', metavar='FILE', help='file containing a list of marker names/identifiers,'
        ' one per line'
    )
    subparser.add_argument(
        '-r', '--region', metavar='RGN', help='restrict results to the '
        'specified genomic region; format chrX:YYYY-ZZZZZ'
    )
    subparser.add_argument(
        '--query', metavar='STRING', help='Retrieve records using a Pandas-style query'
    )
    subparser.add_argument('id', nargs='*', help='marker identifier')


def main(args):
    if args.ae_pop:
        microhapdb.set_ae_population(popid=args.ae_pop)
    if args.GRCh37:
        microhapdb.set_reference(37)
    if args.query:
        if len(args.id) > 0 or args.panel is not None:
            warning = 'WARNING: ignoring user-supplied marker IDs in --query mode'
            print('[MicroHapDB::marker]', warning, file=sys.stderr)
        result = microhapdb.markers.query(args.query, engine='python')
    elif args.region:
        if len(args.id) > 0 or args.panel is not None:
            warning = 'WARNING: ignoring user-supplied marker IDs in --region mode'
            print('[MicroHapDB::marker]', warning, file=sys.stderr)
        result = microhapdb.retrieve.by_region(args.region)
    elif len(args.id) > 0 or args.panel is not None:
        idents = args.id
        if args.panel:
            if hasattr(microhapdb.panel, args.panel):
                func = getattr(microhapdb.panel, args.panel)
                idents.extend(func())
            else:
                with open(args.panel, 'r') as fh:
                    for line in fh:
                        idents.append(line.strip())
        idents = microhapdb.retrieve.standardize_marker_ids(idents)
        result = microhapdb.markers[microhapdb.markers.Name.isin(idents)]
    else:
        result = microhapdb.markers
    if args.format == 'table':
        print_table(result, trunc=args.trunc)
    elif args.format == 'detail':
        print_detail(result, delta=args.delta, minlen=args.min_length, extendmode=args.extend_mode)
    elif args.format == 'fasta':
        print_fasta(result, delta=args.delta, minlen=args.min_length, extendmode=args.extend_mode)
    elif args.format == 'offsets':
        refr = "Hg37" if args.GRCh37 else "Hg38"
        print_offsets(result, delta=args.delta, minlen=args.min_length, extendmode=args.extend_mode, refr=refr)
    else:
        raise ValueError(f'unsupported view format "{args.format}"')
    if args.ae_pop:
        microhapdb.set_ae_population(popid=None)  # Reset
    if args.GRCh37:
        microhapdb.set_reference(38)  # Reset
