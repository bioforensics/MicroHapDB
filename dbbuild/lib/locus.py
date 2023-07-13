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

from collections import defaultdict

class Locus(list):
    def __init__(self):
        self.markers = list()
        self.markers_by_definition = defaultdict(list)
        self.definition_names = dict()
        self.source_name_map = defaultdict(dict)

    def __len__(self):
        return len(self.markers)

    def add(self, marker):
        if len(self) > 0:
            if marker.name != self.markers[-1].name:
                raise ValueError(f"MH name mismatch: {marker.name} vs. {self.markers[-1].name}")
        self.markers.append(marker)
        self.markers_by_definition[marker.posstr()].append(marker)

    @property
    def name(self):
        if len(self) < 1:
            return None
        return self.markers[0].locus

    def resolve(self):
        if len(self) <= 1:
            for marker in self.markers:
                yield marker
            return
        self.check_overlap()
        for marker in sorted(self.markers, key=lambda m: (m.source.year, m.source.order, m.name.lower())):
            if marker.posstr() in self.definition_names:
                message = f"Marker {marker.name} as defined in {marker.sources[0].name} was defined previously and is redundant"
                print(message)
                self.source_name_map[marker.sources[0].name][marker.name] = self.definition_names[marker.posstr()]
                continue
            else:
                new_name = marker.name
                if len(self.markers_by_definition) > 1:
                    new_name = f"{marker.name}.v{len(self.definition_names) + 1}"
                self.definition_names[marker.posstr()] = new_name
                self.source_name_map[marker.sources[0].name][marker.name] = new_name
                marker.name = new_name
                for othermarker in self.markers_by_definition[marker.posstr()]:
                    if othermarker != marker:
                        marker.sources.append(othermarker.sources[0])
                yield marker

    def check_overlap(self):
        if len(self) <= 1:
            return
        for i in range(len(self.markers)):
            marker1 = self.markers[i]
            for j in range(len(self.markers)):
                if i == j:
                    continue
                marker2 = self.markers[j]
                if marker1.overlaps(marker2):
                    break
            else:
                raise ValueError(f"{marker1.name} ({marker1.sourcename}) does not overlap with any other markers defined at {marker1.locus}")
