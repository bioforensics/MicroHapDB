# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

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
