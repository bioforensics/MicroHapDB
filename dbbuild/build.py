#!/usr/bin/env python

from lib import SourceIndex

index = SourceIndex("sources", "databases/dbSNP", "databases")
index.interval_check()
index.update_marker_names()
index.markers.to_csv("marker.csv", index=False)
index.indels.to_csv("indels.csv", index=False)
index.frequencies.to_csv("frequency.csv", index=False, float_format="%.5f")
index.populations.to_csv("population.csv", index=False)
index.merges.to_csv("merged.csv", index=False)
print(index)
