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
from textwrap import dedent


def main(args):
    result = microhapdb.retrieve_by_id(args.id)
    print(result.to_string(index=False))


def subparser(subparsers):
    desc = (
        microhapdb.cli.bubbletext + "\nRetrieve marker or population records by name or identifier"
    )
    epilog = """\
    Examples::

        microhapdb lookup rs10815466
        microhapdb lookup mh12KK-043
        microhapdb lookup Japanese
    """
    epilog = dedent(epilog)
    subparser = subparsers.add_parser(
        "lookup",
        description=desc,
        epilog=epilog,
        formatter_class=RawDescriptionHelpFormatter,
    )
    subparser.add_argument("id", help="record identifier")
