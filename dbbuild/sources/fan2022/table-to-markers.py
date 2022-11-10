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
    with open(args.table, "r") as table, open(args.indels, "w") as fh:
        print("Marker", "VariantIndex", "Refr", "Alt", sep="\t", file=fh)
        for marker, index, refr, alt in parse_indels(table):
            print(marker, index, refr, alt, sep="\t", file=fh)

def get_parser():
    parser = ArgumentParser()
    parser.add_argument("table")
    parser.add_argument("tsv")
    parser.add_argument("bed")
    parser.add_argument("indels")
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


def parse_indels(tsv):
    next(tsv)
    marker = None
    index = -1
    for line in tsv:
        if line.strip() == "":
            continue
        fields = line.rstrip().split("\t")
        if fields[0].startswith("MH"):
            marker = "mh" + fields[0][2:]
            index = -1
        index += 1
        refr = fields[6]
        alts = fields[7].split("/")
        altlengths = sorted([len(alt) for alt in alts])
        if len(refr) > 1 or altlengths[-1] > 1:
            yield marker, index, refr, ",".join(alts)


if __name__ == "__main__":
    main()
