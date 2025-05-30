#!/bin/bash
merged_fastq_dir=/mnt/research/EvansLab/Isabella/merged_fastq
usearch=/mnt/research/EvansLab/Software/usearch_linux_x86_12.0-beta

cd $merged_fastq_dir
$usearch -fastq_filter merged_reads.fastq -fastq_trunclen 380 -fastq_maxee 1.0 -fastaout merged_reads_filtered.fa
