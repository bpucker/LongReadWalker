# LongReadWalker

Inspired by primer walking in the first sequencing projects, this script is intended to close genome assembly gaps by long read walking. The focus is on a short run time and a very small number of depencencies (BLAST).


### Background ###

Initially, tools for gap filling relied on short reads e.g. [Sealer](https://doi.org/10.1186/s12859-015-0663-4), [GapBlaster](https://doi.org/10.1371/journal.pone.0155327), and [GapFiller](https://doi.org/10.1186/1471-2105-13-S14-S8). Now, there are even some tools which can utilize long reads e.g. [GAPPadder](https://doi.org/10.1186/s12864-019-5703-4), [TGS-GapCloser](https://doi.org/10.1093/gigascience/giaa094), [PGcloser](https://doi.org/10.1177%2F1176934320913859), and [gapFinisher](https://doi.org/10.1371/journal.pone.0216885). While some of these tools might not be suitable for the large size of plant genome sequencing datasets, others might require the installation of dependencies. Therefore, I developed the LongReadWalker as a single Python script to close gaps harnessing the power of long sequence reads.


### Concept ###

The LongReadWalker needs a start sequence i.e. a sequence at a contig end. All sequence reads are screened for a match and the read which allows the best extension is used for the first step. The sequence of this read is used to extend the "assembled" sequence and to define a new query sequence at the unmatched read end. This query is searched against all reads for the next extension step. This simple combination of extension and search is repeated for a defined number of iteration. Depending on the read length distribution 10 steps can be sufficient to "assemble" one million basepairs which should be enough to span most gaps.
If both flanking sequences of a gap are known, the gap filling via long read walking can be started from both ends. This allows to compare the results and can be a way to increased the accuracy. The major drawback of this method is the low accuracy of the sequence placed in the gap. Since the gap is closed by sequences taken from the reads, the accuracy depends on the quality of these reads. Due to the relatively high error rate of long reads, the accuracy might be only 85%. Therefore, it is recommend to perform this gap filling prior to the polishing steps, which would also improve the quality of the sequence inserted into the gap.


### Validation and benchmarking ###

LongReadWalker was tested on long read sequencing datasets of different species:
1) Arabidopsis thaliana: Based on a flanking sequence, a fictive gap of 1 Mbp was assembled in 10 steps. This gap filling was performed in both directions and validated through comparison against the A. thaliana TAIR10 reference genome sequence via dot plots.
2) [add upon publication]
3) [add upon publication]


### Reference (how to cite) ###

Please cite this github repository when using LongReadWalker. 
