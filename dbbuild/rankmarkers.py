#!/usr/bin/env python
import microhapdb

m = microhapdb.markers.copy().round({'In': 4})
m['Start'] = m.Offsets.apply(lambda x: min(map(int, x.split(','))))
m['End'] = m.Offsets.apply(lambda x: max(map(int, x.split(','))) + 1)
m['Length'] = m.apply(lambda x: x.End - x.Start, axis=1)

best_by_ae = m\
    .sort_values(['AvgAe'], ascending=False)\
    .reset_index(drop=True)\
    .reset_index()\
    .rename(columns={'index': 'RankAvgAe'})
best_by_ae['RankAvgAe'] = best_by_ae['RankAvgAe'] + 1
best_by_in = m\
    .sort_values(['In'], ascending=False)\
    .reset_index(drop=True)\
    .reset_index()\
    .rename(columns={'index': 'RankIn'})
best_by_in['RankIn'] = best_by_in['RankIn'] + 1

best_by_ae\
    .join(best_by_in[['Name', 'RankIn']].set_index('Name'), on='Name')\
    [['Name', 'AvgAe', 'RankAvgAe', 'In', 'RankIn', 'Length']]\
    .sort_values(['AvgAe'], ascending=False)\
    .head(50)\
    .to_csv('microhaps-highest-ae.tsv', sep='\t', index=False)

best_by_in\
    .join(best_by_ae[['Name', 'RankAvgAe']].set_index('Name'), on='Name')\
    [['Name', 'In', 'RankIn', 'AvgAe', 'RankAvgAe', 'Length']]\
    .sort_values(['In'], ascending=False)\
    .head(50)\
    .to_csv('microhaps-highest-in.tsv', sep='\t', index=False)
