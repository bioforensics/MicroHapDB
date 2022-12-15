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

import pandas as pd


def reformat_markers(infile, outfile):
    final_defs = list()
    markers = pd.read_csv(infile, sep="\t")
    for marker, mdata in markers.groupby("MarkerName"):
        name = "mh" + marker[2:6] + "-" + marker[6:]
        chrom = f"chr{mdata.Chrom.iloc[0]}"
        rsids = mdata.SNPID.tolist()
        offsets = [pos for pos in mdata.Position]
        rsidstr = ";".join(rsids)
        offsetstr = ";".join(map(str, offsets))
        entry = (name, None, len(rsids), "GRCh38", chrom, offsetstr, rsidstr)
        final_defs.append(entry)
    colnames = ["Name", "Xref", "NumVars", "Refr", "Chrom", "Positions", "VarRef"]
    defs = pd.DataFrame(final_defs, columns=colnames)
    defs.to_csv(outfile, index=False)


def reformat_frequencies(infile, outfile):
    freqdata = list()
    freqs = pd.read_csv(infile, sep="\t")
    for n, row in freqs.iterrows():
        standardname = "mh" + row.MarkerName[2:6] + "-" + row.MarkerName[6:]
        haplotype = "|".join(row.Allele)
        entry = (standardname, "MHDBP-48c2cfb2aa", haplotype, row.Frequency)
        freqdata.append(entry)
    freqtable = pd.DataFrame(freqdata, columns=["Marker", "Population", "Allele", "Frequency"])
    freqtable.to_csv(outfile, index=False, float_format="%.4f")
