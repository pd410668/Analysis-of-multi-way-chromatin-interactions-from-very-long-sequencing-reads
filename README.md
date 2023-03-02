# Analysis of multi-way chromatin interactions from very long sequencing reads
The thesis is focused on the concept of “chromosomal walks” introduced by Olivares-Chauvet et
al., 2016. These chromosomal walks are derived from interactions identified de novo in whole-
genome chromatin conformation data. The structure of chromosomal walks is linear, linking
multiple genomic loci together into a chain. In the thesis, a Snakemake workflow to analyze
this type of interactions in a reproducible manner is implemented. Individual scripts are
written in Python. The workflow is applied to published Hi-C data from human K562 cancer
cells and from mouse embryonic stem cells. Selected properties of identified chromosomal
walks, in particular the frequency of staying in the same Topologically Associating Domain,
have been compared.
