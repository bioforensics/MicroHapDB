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


from collections import defaultdict
import pandas as pd


def compile_frequencies(markerdefs, freqs, outcsv):
    with open(markerdefs, "r") as markerfh, open(freqs, "r") as freqfh:
        haplotypes = parse_haplotypes(markerfh)
        freqdata = parse_frequencies(freqfh, haplotypes)
        freqdata.to_csv(outcsv, index=False, float_format="%.03f")


def parse_haplotypes(stream):
    haplotypes = defaultdict(dict)
    marker = None
    for line in stream:
        values = line.strip().split(",")
        if line.startswith("mh"):
            marker = values[0]
        else:
            hapid = values[0][2:]
            allelelist = [a.strip() for a in values[1].split(":")]
            haplotype = "|".join(allelelist)
            haplotypes[marker][hapid] = haplotype
    return haplotypes


def parse_frequencies(stream, haplotypes):
    popids = {
        "Africa": "MHDBP-3dab7bdd14",
        "Asia": "MHDBP-936bc36f79",
        "NL": "MHDBP-383d86606a",
    }
    freqdata = list()
    marker = None
    next(stream)
    for line in stream:
        values = line.strip().split("\t")
        if line.startswith("mh"):
            marker = values[0]
        else:
            if marker == "mh17PK-86511":
                continue
            hapid = values[0]
            freqs = values[1:]
            for freq, pop in zip(freqs, ("NL", "Asia", "Africa")):
                entry = (marker, popids[pop], haplotypes[marker][hapid], float(freq))
                freqdata.append(entry)
    return pd.DataFrame(freqdata, columns=["Marker", "Population", "Allele", "Frequency"])


def compile_marker_definitions(markerrsids, markerdefs, outcsv):
    markers = check_definitions(markerrsids, markerdefs)
    markers.to_csv(outcsv, index=False)


def check_definitions(markerfile, suppfile):
    finaldefs = list()
    marker_rsids = pd.read_csv(markerfile, sep="\t")
    markers = collapse_marker_definitions(marker_rsids)
    with open(suppfile, "r") as fh:
        for name, chrom, astart, aend, offsets in parse_supp_data(fh):
            mdata = markers[markers.Marker == name].iloc[0]
            assert offsets == mdata.Offsets38, mdata
            o38 = ";".join(map(str, mdata.Offsets38))
            rsids = ";".join(mdata.RSIDs)
            entry = (name, None, len(offsets), "GRCh38", chrom, o38, rsids)
            finaldefs.append(entry)
    colnames = ["Name", "Xref", "NumVars", "Refr", "Chrom", "Positions", "VarRef"]
    return pd.DataFrame(finaldefs, columns=colnames)


def collapse_marker_definitions(markers):
    data = list()
    for marker, mdata in markers.groupby("Marker"):
        o38 = [int(v.split(":")[1]) for v in mdata.Position38]
        rs = mdata[mdata.RSID != "-"].RSID.tolist()
        entry = (marker, o38, rs)
        data.append(entry)
    return pd.DataFrame(data, columns=["Marker", "Offsets38", "RSIDs"])


def parse_supp_data(stream):
    for line in stream:
        if not line.startswith("mh"):
            continue
        name, chrom, amplstart, amplend, localoffsets = line.strip().split(",")
        offsets = [int(amplstart) + int(o) + 1 for o in localoffsets.split(":")]
        offsets = sorted(offsets)
        if amplstart > amplend:
            amplstart, amplend = amplend, amplstart
        yield name, chrom, int(amplstart), int(amplend), offsets
