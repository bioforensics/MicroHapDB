# Ranking markers

MicroHapDB provides three primary criteria for ranking markers.
Phased genotypes for 2,504 individuals from Phase 3 of the 1000 Genomes Project, selected from 26 global population groups, are used to calculate these statistics.


## Effective number of alleles (*A<sub>e</sub>*)

The marker's effective number of alleles ($A_e$) is computed individually for 26 population groups.
By default, the average of the 26 populations $A_e$ is shown, but the `--ae-pop` flag or the `microhapdb.set_ae_population` function can be used to specify a single population for which to display $A_e$ values.

The $A_e$ statistic is a measure of the *within-population* allelic variation at a locus, which corresponds to the marker's diagnostic power for identification purposes.


## Informativeness (*I<sub>n</sub>*)

Rosenberg's informativeness for assignment to populations ($I_n$) is computed on 26 global population groups.


## Fixation index (*F<sub>ST</sub>*)

The fixation index ($F_{ST}$) is computed on 26 global population groups.

The $I_n$ and $F_{ST}$ statistics measure *between-population* allelic variation at a locus, which corresponds to the marker's utility for predicting population of origin.
