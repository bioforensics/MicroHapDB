# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from locus import alfred_locus_list


ALL_ALFRED_LOCI = list()
with open('alfred/Microhap_alleleF_198.txt') as instream:
    for locusid, locusname in alfred_locus_list(instream):
        ALL_ALFRED_LOCI.append(locusid)


rule loci:
    input: expand('alfred/downloads/locus-detail/{locus}.html.gz', locus=ALL_ALFRED_LOCI)


rule locusdetail:
    output: 'alfred/downloads/locus-detail/{locus}.html.gz'
    shell: 'curl https://alfred.med.yale.edu/alfred/recordinfo.asp?UNID={wildcards.locus} | gzip -c > {output}'
