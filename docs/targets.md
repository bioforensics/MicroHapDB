# Configuring marker target sequences

When designing an assay for targeted amplification or hybridization enrichment of microhaplotype markers, the precise location of the targeted locus must be specified for each marker.
This typically includes more than just the sequence bounded by the first and last SNPs\* in the marker.
When targeting multiple markers in a single assay, the relative lengths of the target sequences and the position of the SNPs within the target sequences must be considered.

*\*NOTE: Some microhap markers include indel polymorphisms as well as SNPs, but this document will refer only to SNPs for clarity.*

MicroHapDB allows the user to specify parameters for configuring target sequences to meet desired length constraints, and recomputes SNP positions relative to the target sequence for convenience.
Target sequence information can be formatted in several different ways, and can be used not only for assay design but also for bioinformatics analysis of microhap sequence data.

In brief, the target sequence of a marker is determined by including $\delta$ bp of sequence flanking both the marker's 5' and 3' bounding SNPs, and then—if needed—extending further in one or both directions (depending on the extension strategy $E$) to satisfy the minimum length criterion $L$.
Each of these parameters is described in detail below.


## Delta: flanking sequence

```text
    δ                                                                     δ
<=========                                                            =========>
          *                          *                               *
AGCAAACAAAGAACAGTATGTGACAGAGACTGTATCTGGTATGCAAAACTTGAAATACTTACTATCTGCCACTTTACAGA
..........G..........................G...............................C..........
..........G..........................G...............................G..........
..........T..........................C...............................C..........
..........T..........................G...............................C..........
```

The $\delta$ (`--delta`) parameter specifies the number of flanking nucleotides to include in the target sequence beyond the 5' and 3' bounding SNPs.
This flanking sequence is included in the target sequence regardless of the value of the minimum length parameter.
The example above shows the marker mh09USC-9qA, whose 3 SNPs occupy 60 bp from the 5' bounding SNP to the 3' bounding SNP.
Above the marker reference sequence (chr9:70747687-70747747 on GRCh38) asterisk symbols mark the location of the SNPs defined for the marker; below, haplotypes observed in 1000 Genomes Project data are shown.

In this example $\delta=10$, so 10 flanking base pairs are included on either side of the marker's bounding SNPs.


## *L*: minimum length criterion

```text
<-------------------------------------- L=90 -------------------------------------------->
          δ                                                                   δ
<~~~~<=========*                          *                               *=========>~~~~>
AAGAAAGCAAACAAAGAACAGTATGTGACAGAGACTGTATCTGGTATGCAAAACTTGAAATACTTACTATCTGCCACTTTACAGAAAAGT
...............G..........................G...............................C...............
...............G..........................G...............................G...............
...............T..........................C...............................C...............
...............T..........................G...............................C...............
```

The $L$ (`--min-length`) parameter specifies the minimum length for a marker target sequence.
This parameter can provide consistency when designing a panel with markers of different lengths.
If a candidate panel contains both very long (hundreds of bp) and very short (10-20 bp or less) microhap markers, balanced amplification or enrichment is typically easier to achieve with target sequences of similar lengths.
The $L$ parameter can be used to improve uniformity of sequence lengths.

The example above again shows marker mh09USC-9qA.
Assume a minimum length criterion of $L=90$ has been chosen.
Even after $\delta=10$ flanking nucleotides have been included on either side of the bounding SNPs (`=====>`), the reference sequence is still only 80 bp in length: 10 bp short.
The minimum length criterion $L=90$ requires that the sequence be extended an additional five base pairs in each direction (`~~~~>`) to reach 90 bp.


## *E*: extend mode

```text
     <~~~~==========------------------------- E=s ------------------------------==========~~~~>
<~~~~~~~~~========== E=5 -------------------------------------------------------=========>
          <=========------------------------------------------------------- E=3 ==========~~~~~~~~~>
                    *                          *                               *
AAAACAAGAAAGCAAACAAAGAACAGTATGTGACAGAGACTGTATCTGGTATGCAAAACTTGAAATACTTACTATCTGCCACTTTACAGAAAAGTTTGCC
....................G..........................G...............................C....................
....................G..........................G...............................G....................
....................T..........................C...............................C....................
....................T..........................G...............................C....................
```

When extending a target sequence to meet a specified minimum length, MicroHapDB's default behavior is to extend both the 5' end and the 3' end symmetrically.
However, it's possible to extend in only one direction, so that the marker's SNPs are closer to one end than the other.
Depending on the length of target sequences and the assay used, it may be helpful to have the SNPs of interest closer to one side of the target than the other.

Together, the $\delta$ and $E$ parameters parameters can be tuned to configure the placement of the associated SNPs within the target sequence.
When extending marker sequences in one direction, the $\delta$ parameter controls the distance of the SNPs of interest from the near end of the target sequence.

In the example above, marker mh09USC-9qA is shown yet again.
This time, the intervals enclosing three possible target sequences are shown for $L=90$.
The first target sequence corresponds to the default symmetric extend mode $E=s$ (`--extend=mode=symmetric`).
The second target corresponds to 5' only extension $E=5$ (`--extend-mode=5`), while the third target corresponds to 3' only extension $E=3$ (`--extend-mode=3`).

## CLI Example

```
microhapdb marker --panel=markerids.txt --delta=25 --min-length=150 --extend-mode=3 \
    --format=fasta \
    > marker-refr-seqs.fasta
microhapdb marker --panel=markerids.txt --delta=25 --min-length=150 --extend-mode=3 \
    --format=offsets \
    > marker-definitions.tsv
```


## Impacted outputs

MicroHapDB provides several different views for marker information.
This includes the default tabular summary format, a detailed human-readable description of the marker, a view for target sequences in Fasta format, and a table of the position of each marker's SNPs.
The default tabular summary does not include target sequence information, but all other views do.
As a result, it is important to use consistent values of $\delta$, $L$, and $E$ parameters for different views of the same data.

For example, haplotype calling with [MicroHapulator](https://microhapulator.readthedocs.io) requires marker reference sequences in Fasta format and marker definitions (SNP offsets) in tabular format.
These files can be prepared using the `--format=fasta` and `--format=offsets` modes of `microhapdb marker`, respectively.
However, if the two files were not prepared using consistent $\delta$, $L$, and $E$ values, the haplotype calls made with MicroHapulator will be incorrect.
