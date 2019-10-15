# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from argparse import RawDescriptionHelpFormatter
import microhapdb
from textwrap import dedent


def subparser(subparsers):
    desc = 'Retrieve population records by identifier or query'
    epilog = """\
    Examples::

        microhapdb population SA004244O
        microhapdb population Han Japanese Korean
        microhapdb population --format=detail 'Melanesian, Nasioi'
        microhapdb population --query='Source == "10.1016/j.fsigen.2018.05.008"'
    """
    epilog = dedent(epilog)
    subparser = subparsers.add_parser(
        'population', description=desc, epilog=epilog, formatter_class=RawDescriptionHelpFormatter,
    )
    subparser.add_argument('--format', choices=['table', 'detail'], default='table')
    subparser.add_argument(
        '--query', metavar='STRING', help='Retrieve records using a Pandas-style query'
    )
    subparser.add_argument('id', nargs='*', help='population identifier(s)')


def main(args):
    if args.query:
        result = microhapdb.populations.query(args.query)
    elif len(args.id) > 0:
        idents = microhapdb.retrieve.standardize_population_ids(args.id)
        result = microhapdb.populations[microhapdb.populations.ID.isin(idents)]
    else:
        result = microhapdb.populations
    print(result.to_string(index=False))
