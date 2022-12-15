#!/usr/bin/env python

from lib import DataSource
from pathlib import Path

sources = list()
for sourcepath in Path("sources").iterdir():
    if not sourcepath.is_dir():
        continue
    sources.append(DataSource(sourcepath))
sources.sort(key=lambda s: (s.year, s.name.lower()))
print(*sources, sep="\n")
