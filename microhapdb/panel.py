# -----------------------------------------------------------------------------
# Copyright (c) 2019, Battelle National Biodefense Institute.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/microhapdb)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from itertools import combinations
import microhapdb


def panel_alpha():
    '''Initial minimum-effort panel selection.

    There are many factors to consider when designing a panel. This method ignores most of those
    considerations and simply grabs the microhap marker from each autosome with the largest
    effective number of alleles (Ae) averaged across all populations. Only populations from ALFRED
    were considered as these were the most comprehensive published resource at the time.
    '''
    markers = microhapdb.markers.query('Source == "ALFRED"').\
        sort_values('AvgAe', ascending=False).\
        drop_duplicates('Chrom')
    return list(markers.Name)


def panel_alpha_legacy():
    '''The results of the initial alpha panel design'''
    return [
        'mh01KK-117', 'mh10KK-163', 'mh11KK-180', 'mh12CP-008', 'mh13KK-218', 'mh14CP-003',
        'mh15CP-001', 'mh16KK-049', 'mh17CP-001', 'mh18CP-005', 'mh19KK-299', 'mh02KK-134',
        'mh20KK-307', 'mh21KK-320', 'mh22KK-061', 'mh03CP-005', 'mh04KK-030', 'mh05KK-170',
        'mh06KK-008', 'mh07CP-004', 'mh08KK-039', 'mh09KK-157'
    ]


def panel_beta():
    '''First attempt at optimizing panel design.

    Slightly less unsophisticated approach to panel selection than panel_alpha,
    but still ignoring some important considerations. This method focuses on a
    few simple filters and simple operations.
    - discard any microhap not present in ALFRED
    - discard any microhap with an average Ae of less than 2.0
    - discard any microhap that spans more than 250 bp
    - for each chromosome, grab the 3 microhaps with the highest combined
      average Ae such that no 2 microhaps occur within 25 Mb; if this criterion
      is too strict, reduce the distance and then the number of desired
      microhaps from 3 to 2 until a compatible set is selected
    - combine microhaps from all chromosomes and select the top 50 by AvgAe
    '''
    markers = microhapdb.markers.copy()
    markers['Start'] = markers.Offsets.apply(lambda x: min(map(int, x.split(','))))
    markers['End'] = markers.Offsets.apply(lambda x: max(map(int, x.split(','))) + 1)
    markers['Length'] = markers.apply(lambda x: x.End - x.Start, axis=1)
    markerids = set()
    for chromid in markers.Chrom.unique():
        chromloci = markers[
            (markers.Source == 'ALFRED')
            & (markers.Chrom == chromid)
            & (markers.AvgAe > 2.0)
            & (markers.Length <= 250)
        ]

        def trycombos(n=3, dist=25e6):
            opt_ae, opt_markers = None, None
            for testmarkerids in combinations(chromloci.Name, n):
                testmarkers = markers[markers.Name.isin(testmarkerids)]
                for coord1, coord2 in combinations(testmarkers.Start, 2):
                    if abs(coord1 - coord2) < dist:
                        break
                else:
                    ae = sum(testmarkers.AvgAe) / len(testmarkers.AvgAe)
                    if opt_ae is None or ae > opt_ae:
                        opt_ae = ae
                        opt_markers = testmarkerids
            return opt_markers
        params = (
            (3, 25e6), (3, 20e6), (2, 25e6), (2, 20e6),
            (3, 15e6), (2, 15e6),
            (3, 10e6), (2, 10e6), (2, 7.5e6),
        )
        for n, dist in params:
            testmarkers = trycombos(n=n, dist=dist)
            if testmarkers is not None:
                markerids.update(testmarkers)
                break
    panel = markers[markers.Name.isin(markerids)].sort_values('AvgAe').head(50)
    return list(panel.Name)


def panel_beta_legacy():
    '''The results of the initial beta panel design'''
    return [
        'mh01CP-016', 'mh01KK-117', 'mh01KK-205', 'mh02KK-136', 'mh02KK-138', 'mh03CP-005',
        'mh03KK-007', 'mh03KK-150', 'mh04KK-013', 'mh04KK-017', 'mh04KK-030', 'mh05KK-020',
        'mh06CP-007', 'mh06KK-008', 'mh07CP-004', 'mh07KK-030', 'mh07KK-031', 'mh08KK-032',
        'mh08KK-039', 'mh09KK-033', 'mh09KK-153', 'mh09KK-157', 'mh10CP-003', 'mh10KK-101',
        'mh10KK-170', 'mh11KK-037', 'mh11KK-180', 'mh11KK-191', 'mh12CP-008', 'mh12KK-046',
        'mh12KK-202', 'mh13KK-213', 'mh13KK-223', 'mh14CP-003', 'mh14KK-048', 'mh15KK-069',
        'mh15KK-104', 'mh16KK-255', 'mh17CP-001', 'mh17KK-052', 'mh17KK-055', 'mh18CP-005',
        'mh18KK-293', 'mh19KK-057', 'mh19KK-299', 'mh20KK-058', 'mh20KK-307', 'mh21KK-315',
        'mh22KK-060', 'mh22KK-061'
    ]
