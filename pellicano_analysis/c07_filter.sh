#!/bin/bash
merged_fastq_dir=/mnt/research/EvansLab/Isabella/nifH_tools/pellicano_analysis/mock_c07
usearch=/mnt/research/EvansLab/Software/usearch_linux_x86_12.0-beta

cd $merged_fastq_dir
$usearch -fastq_filter c07_merged_fastq.fq -fastq_trunclen 380 -fastq_maxee 1.0 -fastaout c07_merged_reads_filtered.fa