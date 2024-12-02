# -------------------------------------------------------------------------------------------------
# Copyright (c) 2024, DHS.
#
# This file is part of MicroHapDB (http://github.com/bioforensics/MicroHapDB) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------


class FilteredMarker:
    field_names = [
        "Marker",
        "Extent",
        "Ae",
        "FailMode",
        "MaxAltFreq",
        "LowComplexRepeats",
        "Repeats",
        "STRs",
    ]

    def __init__(self, mask_data):
        self.mask_data = mask_data
        for key, entry in mask_data.items():
            if entry is None:
                continue
            self.marker = entry.name
            self.extent = entry.Extent
            self.ae = entry.Ae
            break
        else:
            raise ValueError("empty mask data?")

    @property
    def fields(self):
        return (
            self.marker,
            self.extent,
            self.ae,
            self.fail_mode,
            self.max_alt_freq,
            self.low_complex_repeats,
            self.repeats,
            self.strs,
        )

    @property
    def fail_mode(self):
        failures = set()
        for key, mask in self.mask_data.items():
            if mask is not None:
                failures.add(key)
        return "".join(f"[{key}]" for key in sorted(failures))

    @property
    def max_alt_freq(self):
        if self.mask_data["indel"] is None:
            return ""
        return self.mask_data["indel"].MaxAltFreq

    @property
    def low_complex_repeats(self):
        if self.mask_data["lowcomplex"] is None:
            return ""
        return self.mask_data["lowcomplex"].Repeats

    @property
    def repeats(self):
        if self.mask_data["repeat"] is None:
            return ""
        return self.mask_data["repeat"].Repeats

    @property
    def strs(self):
        if self.mask_data["str"] is None:
            return ""
        return self.mask_data["str"].Locus

    @property
    def filtering(self):
        return (
            self.mask_data["indel"] is not None,
            self.mask_data["length"] is not None,
            self.mask_data["lowcomplex"] is not None,
            self.mask_data["repeat"] is not None,
            self.mask_data["str"] is not None,
        )
