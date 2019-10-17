# -----------------------------------------------------------------------------
# Copyright (c) 2018, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from marker import parse_markers_from_table


MARKERS = parse_markers_from_table()


rule markers:
    input: expand('downloads/marker-detail/{marker}.html.gz', marker=MARKERS)


rule frequencies:
    output: 'downloads/Microhap_alleleF_198.txt'
    shell: 'curl https://alfred.med.yale.edu/alfred/selectDownload/Microhap_alleleF_198.txt > {output}'


rule markerdetail:
    output: 'downloads/marker-detail/{marker}.html.gz'
    shell: 'curl https://alfred.med.yale.edu/alfred/recordinfo.asp?UNID={wildcards.marker} | gzip -c > {output}'
