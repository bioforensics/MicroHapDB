#!/usr/bin/env python

from argparse import ArgumentParser
from collections import namedtuple


Marker = namedtuple("Marker", "name nvars chrom off38 varref")


def main():
    args = get_parser().parse_args()
    header = ["Name", "Xref", "NumVars", "Chrom", "OffsetsHg37", "OffsetsHg38", "VarRef"]
    with open(args.table, "r") as table, open(args.tsv, "w") as tsv, open(args.bed, "w") as bed:
        print(*header, sep="\t", file=tsv)
        for marker in parse_markers(table):
            print_marker_tsv(marker, tsv)
            print_marker_bed(marker, bed)


def get_parser():
    parser = ArgumentParser()
    parser.add_argument("table")
    parser.add_argument("tsv")
    parser.add_argument("bed")
    return parser


def parse_markers(tsv):
    next(tsv)
    marker = None
    for line in tsv:
        if line.strip() == "":
            continue
        fields = line.rstrip().split("\t")
        off38 = int(fields[4]) - 1
        varref = fields[8] if len(fields) >= 9 else ""
        if fields[0].startswith("MH"):
            if marker is not None:
                yield marker
            name = "mh" + fields[0][2:]
            nvars = int(fields[1])
            chrom = int(fields[2][3:])
            varreflist = list() if varref == "" else [varref]
            marker = Marker(name, nvars, chrom, [off38], varreflist)
        else:
            marker.off38.append(off38)
            if varref != "":
                marker.varref.append(varref)
    yield marker


def print_marker_tsv(marker, fh):
    print(
        marker.name,
        "",
        marker.nvars,
        marker.chrom,
        "",
        ",".join(map(str, marker.off38)),
        ",".join(marker.varref),
        sep="\t",
        file=fh,
    )


def print_marker_bed(marker, fh):
    for offset in marker.off38:
        print(f"chr{marker.chrom}", offset, offset + 1, file=fh)


if __name__ == "__main__":
    main()
