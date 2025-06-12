#!/bin/bash
merged_fastq_dir=/mnt/research/EvansLab/Isabella/nifH_tools/pellicano_analysis/mock_c07
usearch=/mnt/research/EvansLab/Software/usearch_linux_x86_12.0-beta

cd $merged_fastq_dir
$usearch -fastq_filter c07_merged_fastq.fq -fastq_trunclen 380 -fastq_maxee 1.0 -fastaout c07_merged_reads_filtered.fa


#results#
00:27 2.1Mb 100.0% Filtering 31.9k passed (89.4%)
00:27 2.1Mb ......      35695  Reads (35.7k)
00:27 2.1Mb ......         48  Discarded reads length < 380
00:27 2.1Mb ......       3719  Discarded reads with expected errs > 1.00
00:27 2.1Mb ......      31928  Filtered reads (31.9k, 89.4%)
00:27 2.1Mb ...... Finished
