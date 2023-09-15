# -------------------------------------------------------------------------------------------------
# Copyright (c) 2023, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

import json
import os
import pandas as pd
from tqdm import tqdm
from urllib.request import urlretrieve

accessions = {37: "GCF_000001405.25", 38: "GCF_000001405.40"}


def progress(t):
    """Stolen shamelessly from the tqdm documentation"""
    last_b = [0]

    def inner(b=1, bsize=1, tsize=None):
        if tsize is not None:
            t.total = tsize
        t.update((b - last_b[0]) * bsize)
        last_b[0] = b

    return inner


rule all:
    input:
        expand(
            "dbSNP/dbSNP_GRCh{version}.{extension}",
            version=(37, 38),
            extension=("vcf.gz", "vcf.gz.tbi"),
        ),
        "dbSNP/refsnp-merged.csv.gz",
        expand("{contrast}.over.chain.gz", contrast=("hg19ToHg38", "hg38ToHg19")),


rule dbsnp:
    output:
        path="dbSNP/dbSNP_GRCh{version}.{extension}",
    wildcard_constraints:
        version="\d+",
    run:
        accession = accessions[int(wildcards.version)]
        url_ext = wildcards.extension.replace("vcf.", "")
        url = f"https://ftp.ncbi.nih.gov/snp/archive/b156/VCF/{accession}.{url_ext}"
        origpath = Path("dbSNP") / Path(url).name
        with tqdm(unit="B", unit_scale=True, leave=True, miniters=1, desc=str(origpath)) as t:
            urlretrieve(url, origpath, reporthook=progress(t))
        os.symlink(origpath.name, Path(output.path).name, dir_fd=os.open("dbSNP", os.O_RDONLY))


rule merged:
    output:
        path="dbSNP/refsnp-merged.csv.gz",
    run:
        url = "https://ftp.ncbi.nih.gov/snp/archive/b156/JSON/refsnp-merged.json.bz2"
        json_path = "dbSNP/refsnp-merged.json.bz2"
        urlretrieve(url, json_path)
        shell(f"bunzip2 {json_path}")
        merged_rsids = dict()
        updateint = 1e6
        threshold = updateint
        with open("dbSNP/refsnp-merged.json", "r") as instream:
            for n, line in enumerate(instream):
                try:
                    data = json.loads(line)
                except Exception:
                    warn(f"Could not parse line {n+1}, skipping: {line}")
                    continue
                source = data["refsnp_id"]
                targets = data["merged_snapshot_data"]["merged_into"]
                for target in targets:
                    merged_rsids[f"rs{source}"] = f"rs{target}"
                if n >= threshold:
                    threshold += updateint
                    if threshold == updateint * 10:
                        updateint = threshold
                    print(f"processed {n} rows")
        table = pd.DataFrame(merged_rsids.items(), columns=["Source", "Target"])
        table.to_csv(output.path, index=False, compression="gzip")


rule chain:
    output:
        path="{contrast}.over.chain.gz",
    run:
        source = wildcards.contrast[:4]
        url = f"https://hgdownload.cse.ucsc.edu/goldenpath/{source}/liftOver/{wildcards.contrast}.over.chain.gz"
        urlretrieve(url, output.path)
