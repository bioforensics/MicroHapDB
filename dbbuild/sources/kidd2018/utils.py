# -------------------------------------------------------------------------------------------------
# Copyright (c) 2018, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

import builtins
from collections import defaultdict
from gzip import open as gzopen
import pandas as pd
from pathlib import Path
from re import findall, search
import sys


def parse_markers_from_table(filename):
    if not Path(filename).is_file():
        return list()
    with smartopen(filename, "r") as fh:
        return [label for label, xref in alfred_marker_list(fh)]


def smartopen(filename, mode):
    """Smart file handler

    Determines whether to create a compressed or un-compressed file handle
    based on filename extension.
    """
    if mode not in ("r", "w"):
        raise ValueError('invalid mode "{}"'.format(mode))
    if filename in ["-", None]:
        filehandle = sys.stdin if mode == "r" else sys.stdout
        return filehandle
    openfunc = builtins.open
    if filename.endswith(".gz"):
        openfunc = gzopen
        mode += "t"
    return openfunc(filename, mode)


def compile_marker_coords(csv, vcf):
    data = list()
    markers = pd.read_csv(csv)
    rsidcoords = dbsnp_to_rsid_coords(vcf)
    for n, row in markers.iterrows():
        name = row.Name
        if name == "mh05KK-058":
            name = "mh15KK-058"
        rsids = row.rsIDs.split(";")
        chrom = rsidcoords[rsids[0]][0]
        positions = [rsidcoords[rsid][1] for rsid in rsids]
        assert len(rsids) == len(positions), name
        posstr = ";".join(map(str, positions))
        data.append((name, row.Xref, len(positions), "GRCh38", chrom, posstr, row.rsIDs))
    colnames = ["Name", "Xref", "NumVars", "Refr", "Chrom", "Positions", "VarRef"]
    return pd.DataFrame(data, columns=colnames)


def dbsnp_to_rsid_coords(vcf):
    rsidcoords = dict()
    with smartopen(vcf, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrnum, posstr, rsid, *values = line.strip().split()
            chrom = chrnum if chrnum.startswith("chr") else "chr" + chrnum
            pos = int(posstr)
            rsidcoords[rsid] = (chrom, pos)
    return rsidcoords


def compile_marker_rsids(filelist):
    data = list()
    for filename in filelist:
        with smartopen(filename, "r") as fh:
            name, rsids = alfred_marker_detail_scrape(fh)
            xref = filename.split("/")[-1].split(".")[0]
            rsidstr = ";".join(rsids)
            data.append((name, xref, rsidstr))
    return pd.DataFrame(data, columns=["Name", "Xref", "rsIDs"])


def alfred_marker_list(instream):
    data = instream.read()
    matches = findall(r"(\S+) \| null \| \d+ \| (\S+)", data)
    for label, xref in matches:
        yield label, xref


def alfred_marker_detail_scrape(instream):
    """Scrape a marker detail page for variant info"""
    data = instream.read()
    vmatch = search(r"in this (\d+)-site microhap", data)
    assert vmatch, instream.name
    nvariants = int(vmatch.group(1))
    pattern = (
        r"(rs\d+)\s*\((A|C|G|T|Ins|Del|Indel -)/(A|C|G|T|Ins|Del|Indel -)"
        r"( *SNP)*\)[\s\S]*?\(UID= *(\S+)\)"
    )
    matches = findall(pattern, data)
    if len(matches) != nvariants:
        message = "mismatch: expected {:d} variants".format(nvariants)
        message += ", only found {:d}".format(len(matches))
        raise ValueError(message)
    dbsnpids = [m[0] for m in matches]
    match = search("(mh\d+(KK|CP|NK)-\d+)", data)
    assert match, data
    name = match.group(1)
    return name, dbsnpids


def compile_pop_data(table, xrefs):
    mapping = pd.read_csv(xrefs, sep="\t")
    data = list()
    with smartopen(table, "r") as fh:
        for popid, name, xref in alfred_pop_data(fh, mapping):
            data.append((popid, name, xref))
    return pd.DataFrame(data, columns=["ID", "Name", "Xref"])


def alfred_pop_data(instream, mapping):
    pops_1kgp = set(mapping.ALFRED)
    popdata = dict()
    for line in instream:
        if line.startswith(("----------", "SI664", "popName")):
            continue
        values = line.strip().split("\t")
        popmatch = search(r"^([^\(]+)\((\S+)\)", values[0])
        assert popmatch, values[0]
        popname = popmatch.group(1)
        label = popmatch.group(2)
        if label in pops_1kgp:
            continue
        typed_sample_size = int(values[1])
        if label in popdata:
            assert popname == popdata[label]
        else:
            popdata[label] = popname
            yield label, popname, ""


def parse_popid_mapping(mappingfile):
    popdict = dict()
    mapping = pd.read_csv(mappingfile, sep="\t")
    for n, row in mapping.iterrows():
        popdict[row.ALFRED] = row.ID1KGP
    return popdict


def alfred_frequencies(table, mappingfile):
    data = list()
    mapping = parse_popid_mapping(mappingfile)
    with smartopen(table, "r") as fh:
        for marker, pop, allele, freq in parse_freqs(fh, mapping):
            if marker == "mh05KK-058":
                marker = "mh15KK-058"
            data.append((marker, pop, allele, freq))
    return pd.DataFrame(data, columns=["Marker", "Population", "Allele", "Frequency"])


def parse_freqs(stream, mapping):
    indels = {
        "SI664579L": {"D": "T", "I": "TA"},
        "SI664597L": {"D": "T", "I": "TG"},
        "SI664640A": {"D": "A", "I": "AATAATT"},
    }

    def cleanup(allelestr, markerid, indels):
        allelestr = allelestr.replace("-", ",")
        if "D" in allelestr or "I" in allelestr:
            allelestr = allelestr.replace("D", indels[markerid]["D"])
            allelestr = allelestr.replace("I", indels[markerid]["I"])
        return allelestr

    line = next(stream)
    if line.startswith("Created on"):
        next(stream)
    for chunk in stream.read().split("-----------------\n"):
        lines = iter(chunk.split("\n"))
        header1 = next(lines)
        xref, _, _, markerid = header1.split(" | ")
        if "-" not in markerid:
            markerid = markerid[:6] + "-" + markerid[6:]
        header2 = next(lines)
        alleles = header2.split()[2:]
        alleles = [cleanup(a, xref, indels) for a in alleles]
        for line in lines:
            if line.strip() == "":
                continue
            values = line.split("\t")
            popid = search(r"^([^\(]+)\((\S+)\)", values[0]).group(2)
            freqs = values[2:]
            if len(alleles) != len(freqs):
                message = "WARNING: allele/freq mismatch "
                message += "for locus " + markerid
                message += " and population " + popid
                message += f"; {len(freqs)} frequencies vs {len(alleles)} alleles"
                raise ValueError(message)
            for allele, freq in zip(alleles, freqs):
                allele = allele.replace(",", "|")
                if popid in mapping:
                    popid = mapping[popid]
                yield markerid, popid, allele, float(freq)
