# ONT-cappable-seq2

We describe here a general overview of the updated ONT-cappable-seq data analysis pipeline.


Data generation: library preparation and sequencing
Total RNA was enriched for primary transcripts using an adapted version of the (SMRT)-cappable-seq enrichment protocol. In parallel, control samples were not enriched but subjected to similar environmental conditions. The RNA was reverse transcribed, PCR amplified and barcoded according to Oxford Nanopore Technology cDNA-PCR protocol (SQK-PCS109 combined with SQK-PBK004). The amplified cDNA samples were pooled together in a final library, which was subsequently loaded on a promethION flowcell (R9.4.1) and sequenced with live base-calling and demultiplexing enabled. Fastq files with Q-scores >7 were retained for further analysis. The overall quality of the sequencing run was assessed using NanoComp (v1.11.2).

Data processing
Note: many of the steps outlined in this workflow are based on the microbepore pipeline that was used for the analysis of prokaryotic Nanopore RNA-seq data (https://felixgrunberger.github.io/microbepore/), with some modifications tailored to our ONT-cappable-seq approach.

I. Data navigation
