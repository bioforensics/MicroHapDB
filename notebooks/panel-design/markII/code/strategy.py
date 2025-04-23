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

from collections import Counter
from filtered import FilteredMarker
import matplotlib
from matplotlib import pyplot as plt
import pandas as pd
from upsetplot import UpSet


class MaskingStrategy(dict):
    def __init__(self, indel_mask, lowcomplex_mask, repeat_mask, str_mask, max_length=260):
        self.masks = {
            "indel": indel_mask,
            "length": None,
            "lowcomplex": lowcomplex_mask,
            "repeat": repeat_mask,
            "str": str_mask,
        }
        self.max_length = max_length

    @classmethod
    def from_tsv(cls, indel_path, lowcomplex_path, repeat_path, str_path, max_length=260):
        return cls(
            pd.read_csv(indel_path, sep="\t").set_index("Marker"),
            pd.read_csv(lowcomplex_path, sep="\t").set_index("Marker"),
            pd.read_csv(repeat_path, sep="\t").set_index("Marker"),
            pd.read_csv(str_path, sep="\t").set_index("Marker"),
            max_length=max_length,
        )

    def __iter__(self):
        for marker in sorted(self.filtered_markers):
            yield self.annotate(marker)

    def apply_masks(self, markers):
        self.masks["length"] = markers[markers.Extent > self.max_length].set_index("Name")
        return markers[~markers.Name.isin(self.filtered_markers)]

    @property
    def filtered_markers(self):
        filtered = set()
        for mask in self.masks.values():
            filtered.update(mask.index)
        return filtered

    @property
    def masked_markers(self):
        data = [self.annotate(mh).fields for mh in sorted(self.filtered_markers)]
        return pd.DataFrame(data, columns=FilteredMarker.field_names)

    def annotate(self, marker_name):
        results = {}
        for key, mask in self.masks.items():
            results[key] = mask.loc[marker_name] if marker_name in mask.index else None
        return FilteredMarker(results)

    def visualize_masking(self, output_path):
        fail_modes = Counter([mh.filtering for mh in self])
        fail_tuples = sorted(fail_modes.items())
        keys = (k for k, v in fail_tuples)
        counts = (v for k, v in fail_tuples)
        index = pd.MultiIndex.from_tuples(
            keys,
            names=[
                "Indel Interference",
                "Length",
                "Low Complexity",
                "Repetitive Content",
                "STR Proximity",
            ],
        )
        mask_data = pd.Series(counts, index=index)
        backend = matplotlib.get_backend()
        plt.switch_backend("Agg")
        figure = plt.figure(figsize=(12, 6), dpi=300)
        UpSet(mask_data, show_counts=True).plot(fig=figure)
        plt.savefig(output_path, bbox_inches="tight")
        plt.switch_backend(backend)
