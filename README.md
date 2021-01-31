# LongReadWalker

Inspired by primer walking in the first sequencing projects, this script is intended to close genome assembly gaps by long read walking. The focus is on a short run time and a very small number of depencencies (BLAST).


### Background ###

Initially, tools for gap filling relied on short reads e.g. [Sealer](https://doi.org/10.1186/s12859-015-0663-4) and GapFiller(https://doi.org/10.1186/1471-2105-13-S14-S8). Now, there are even some tools which can utilize long reads e.g. [GAPPadder](https://doi.org/10.1186/s12864-019-5703-4), [TGS-GapCloser](https://doi.org/10.1093/gigascience/giaa094), [PGcloser](https://doi.org/10.1177%2F1176934320913859), and [gapFinisher](https://doi.org/10.1371/journal.pone.0216885). While some of these tools might not be suitable for the large size of plant genome sequencing datasets, others might require the installation of dependencies. Therefore, I developed the LongReadWalker as a single Python script to close gaps harnessing the power of long sequence reads.

### Concept ###

The LongReadWalker needs a start sequence i.e. a sequence at a contig end. All sequence reads are screened for a match and the read which allos the best extension is used for the first step. The unmatched end of this read is searched against all reads for the next extension step. This simple combination of extension and search is repeated for a defined number of iteration. Depending on the read length distribution 10 steps can be sufficient to "assemble" one million basepairs which should be enough to span most gaps.


