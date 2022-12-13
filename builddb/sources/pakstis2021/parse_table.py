#!/usr/bin/env python

from argparse import ArgumentParser
from collections import namedtuple
import pandas as pd


MarkerInfo = namedtuple("MarkerInfo", "name chrom nsnps off37 rsids")


def main(args):
    table = list()
    with open(args.intable, "r") as fh:
        for marker in parse_markers(fh):
            offsets = ";".join(map(str, marker.off37))
            rsids = ";".join(marker.rsids)
            table.append((marker.name, None, marker.nsnps, "GRCh37", marker.chrom, offsets, rsids))
    colnames = ["Name", "Xref", "NumVars", "Refr", "Chrom", "Positions", "VarRef"]
    data = pd.DataFrame(table, columns=colnames)
    data.to_csv(args.outtable, index=False)


def parse_markers(infile):
    marker = None
    next(infile)
    for line in infile:
        line = line.rstrip()
        if line == "":
            continue
        fields = line.split()
        if len(fields) == 8:
            if marker is not None:
                assert marker.nsnps == len(marker.off37)
                yield marker
            offset = int(fields[4]) - 1
            rsid = fields[5]
            marker = MarkerInfo(
                name=fields[0],
                chrom=f"chr{fields[1]}",
                nsnps=int(fields[3]),
                off37=[offset],
                rsids=[],
            )
            if rsid != ".":
                marker.rsids.append(rsid)
        else:
            offset = int(fields[0]) - 1
            rsid = fields[1]
            marker.off37.append(offset)
            if rsid != ".":
                marker.rsids.append(rsid)
    assert marker.nsnps == len(marker.off37), marker.name
    yield marker


def get_parser():
    parser = ArgumentParser()
    parser.add_argument("intable")
    parser.add_argument("outtable")
    return parser


if __name__ == "__main__":
    main(get_parser().parse_args())
