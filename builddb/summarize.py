#!/usr/bin/env python

from lib import SourceIndex

index = SourceIndex("sources", "databases/dbSNP", "databases")
index.update_marker_names()
markers = index.marker_definitions()
markers.to_csv("marker.csv", index=False)
print(index)
