#!/usr/bin/env python

from argparse import ArgumentParser
from collections import defaultdict
import pandas as pd


def main():
    args = get_parser().parse_args()
    mapping = load_mapping(args.bed37, args.bed38)
    table = pd.read_csv(args.tsv, sep="\t")
    update_table(table, mapping)
    table.to_csv(args.marker, sep="\t", index=False)


def get_parser():
    parser = ArgumentParser()
    parser.add_argument("tsv")
    parser.add_argument("bed37")
    parser.add_argument("bed38")
    parser.add_argument("marker")
    return parser


def load_mapping(bed37, bed38):
    mapping = defaultdict(dict)
    with open(bed37, "r") as fh37, open(bed38, "r") as fh38:
        for line37, line38 in zip(fh37, fh38):
            chromstr37, pos37, _ = line37.split()
            chromstr38, pos38, _ = line38.split()
            chrom37, chrom38 = chromstr37[3:], chromstr38[3:]
            assert chrom38 == chrom37
            mapping[int(chrom38)][int(pos38)] = int(pos37)
    return mapping


def update_table(table, mapping):
    for i, row in table.iterrows():
        off38 = row.OffsetsHg38
        off38 = [int(o) for o in off38.split(",")]
        off37 = [mapping[row.Chrom][o] for o in off38]
        off37 = ",".join(map(str, sorted(off37)))
        assert off37.count(",") + 1 == row.NumVars, (row.Name, row.NumVars, off37.count(","))
        table.loc[i, "OffsetsHg37"] = off37
        table.loc[i, "Chrom"] = f"chr{row.Chrom}"


if __name__ == "__main__":
    main()
