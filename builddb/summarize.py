#!/usr/bin/env python

from lib import SourceIndex

index = SourceIndex("sources", "databases/dbSNP", "databases", exclude=["Auton2015"])
index.update_marker_names()
index.markers.to_csv("marker.csv", index=False)
index.indels.to_csv("indels.csv", index=False)
index.frequencies.to_csv("frequency.csv", index=False, float_format="%.5f")
print(index)
