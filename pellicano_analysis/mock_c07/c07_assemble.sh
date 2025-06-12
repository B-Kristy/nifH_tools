#!/bin/bash

usearch=/mnt/research/EvansLab/Software/usearch_linux_x86_12.0-beta
usearch -fastq_mergepairs Nif_mock_PL4_C07_S338_L001_R1_001.fastq -reverse Nif_mock_PL4_C07_S338_L001_R2_001.fastq -fastqout c07_merged_fastq.fq

# Count the number of total reads before filtering 
wc -l merged_reads.fastq | awk '{print $1/4}' # 35,695 reads