#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from gzip import open as gzopen
import pandas as pd


def collate_marker_definitions(markerrsids, vcf, outcsv):
    with gzopen(vcf, "rt") as fh:
        rsidcoords = parse_rsid_coords(fh)
    markers = marker_coords(markerrsids, rsidcoords)
    markers.to_csv(outcsv, index=False)


def parse_rsid_coords(vcf):
    rsidcoords = dict()
    for line in vcf:
        if line.startswith("#"):
            continue
        chrnum, posstr, rsid, *values = line.strip().split()
        chrom = chrnum if chrnum.startswith("chr") else "chr" + chrnum
        pos = int(posstr)
        rsidcoords[rsid] = (chrom, pos)
    return rsidcoords


def marker_coords(markerfile, rsidcoords):
    markers = pd.read_csv(markerfile, sep="\t")
    markers = markers[markers.Marker != "mh04NH-08"]  # Problematic marker definition; see README
    table = list()
    for n, row in markers.iterrows():
        rsids = row.Variants.split(",")
        offsets = [rsidcoords[r][1] for r in rsids]
        offsetstr = ";".join(map(str, offsets))
        variants = row.Variants.replace(",", ";")
        marker = (row.Marker, None, len(offsets), "GRCh38", row.Chrom, offsetstr, variants)
        table.append(marker)
    colnames = ["Name", "Xref", "NumVars", "Refr", "Chrom", "Positions", "VarRef"]
    return pd.DataFrame(table, columns=colnames)
