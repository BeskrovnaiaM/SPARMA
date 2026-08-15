# SPARMA

This repository contains a script for analyzing the motifs of double-strand breaks created by the nuclease activity of the SPARMA defense system.

## Installation

You can clone this repository or download the source code. 

##### Requirements:

python=3.12.3

bedtools=2.31.1

weblogo=3.7.9

pysam=0.23.3


## `analyze_DSB_motif.py`

This script identifies double-strand break (DSB) sites from + and - reads in a BAM file and analyzes the sequence motif surrounding them. It extracts the 5′-end coordinates of primary alignments separately for + and - strand reads, then uses binary search to find all read pairs whose 5′-ends fall within a specified window. The distance is calculated between 5'-ends of reads on the + and - strands. Positive values indicate that the reads overlap, while negative values indicate that they are separated by the specified number of nucleotides. Reads whose pair distance equals as defined target are flagged and used to build a BED file of their surrounding genomic windows. These regions are extracted as FASTA sequences with bedtools from reference FASTA file and passed to WebLogo to generate a sequence motif logo centered on the break site. The script also tabulates the full distribution of pairwise distances between + and - strand reads and writes it out as a separate file. 
