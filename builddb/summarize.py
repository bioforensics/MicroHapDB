#!/usr/bin/env python

from lib import SourceIndex
from pathlib import Path

index = SourceIndex("sources", "databases/dbSNP")
# print(index)
markers = index.marker_definitions()
markers.to_csv("marker.csv", index=False)
