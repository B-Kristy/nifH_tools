#!/bin/bash

module purge
module load Miniforge3

merged_fastq_dir=/mnt/research/EvansLab/Isabella/merged_fastq/
usearch=/mnt/research/EvansLab/Software/usearch_linux_x86_12.0-beta


# Make an OTU table based on filter corrected seqs
$usearch -otutab $merged_fastq_dir/merged_reads_filtered_hmm.fa -otus $merged_fastq_dir/nifH_corr_nucl_only_nifH.fasta -otutabout otu_table.txt -mapout map.txt -threads 32
