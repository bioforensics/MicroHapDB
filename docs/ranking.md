# Ranking markers

The literature on microhaplotypes describe several measures of allelic variation with which to compare microhap markers.
These include, among others, the following statistics.

- Effective number of alleles ($A_e$)
- Informativeness for assignment ($I_n$)
- Fixation index ($F_{ST}$)

Simply put, $A_e$ measures *within-population* allelic variation, while $I_n$ and $F_{ST}$ measure *between-population* allelic variation.

In the past, MicroHapDB pre-computed all three statistics for each marker in the database.
However, because the values are correlated, and to better handle the increasing scale of the database, MicroHapDB now reports only the $A_e$ for each marker.
$A_e$ calculations are based on the most recent 1000 Genomes Project haplotypes and are provided for 26 global population groups, 4 (non-admixed) continental superpopulation groups, and a global average.
