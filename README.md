[![DOI](https://zenodo.org/badge/332404739.svg)](https://zenodo.org/badge/latestdoi/332404739)

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


### Usage ###


```
Usage:
  python LRW.py --reads <FILE> -seed <FILE> --out <DIR>

Mandatory:
  --reads     STR       FASTA/FASTQ file containing reads
  --seed      STR       FASTA file containing start sequence
  --out       STR       Output directory

Optional:
  --rounds    INT       Rounds of extension [5]
  --block     INT       Query block size [2000]
  --direction STR       Direction of extension (up|down) [down]
  --simcut    INT       Minimal BLAST hit similarity [80]
  --lencut    INT       Minimal BLAST hit length [500]
  --threads   INT       Number of threads for BLAST [8]
```
				


`--reads` specifies an input file in FASTA or FASTQ format containing the reads for the gap filling. gzip compressed files are supported and recognized by the file extension (.gz or .gzip). The file type is checked by reading the first few lines looking for '>' and '@' at the starts.

`--seed` specifies an input file in FASTA format which contains the starting sequence (seed). This sequence should be the end of a contig where the extesion via long read walking should start.

`--out` specifies an output folder where all result files and temporary files will be stored. This folder will be created if it does not exist already.

`--rounds` specifies the number of extension rounds to perform. The number of required rounds can be estimated based on the read length distribution and the expected gap size. When ONT long reads are used, 10 rounds might be sufficient to cover 1 Mbp. Default value: 5.

`--block` species the length of the query block used to find the next read (unit is basepairs). It is important to note that `--lencut` must always be smaller than `--block` to find overlaps. Default: 2000.

`--direction` specifies if the gap is located upstream or downstream of the given sequence. Default is extension towards downstream (down).

`--simcut` specifies the minimal similarity of BLAST hits to be considered for the extension. This value must be smaller than the accuracy of the reads. Latest improvements in basecalling software might allow to increase this value up to 90%. Default: 80%.

`--lencut` specifies the minimal alignment length of BLAST hits to be considered for the extension. This parameter should not be used for stringent filtering, because BLAST was developed to find local alignments. Therefore, it is likely that the best match will not cover the entire query in one block. Default: 500bp.

`--threads` specifies the number of threads used for the BLAST search. This step determines the run time and is repeated in each round. Default: 8 CPUs.




### Reference (how to cite) ###

Please cite this github repository when using LongReadWalker. 
