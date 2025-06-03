#!/bin/bash
# First, dereplicate the sequences
module purge
module load Miniforge3

merged_fastq_dir=/mnt/research/EvansLab/Isabella/merged_fastq
usearch=/mnt/research/EvansLab/Software/usearch_linux_x86_12.0-beta

# Dereplicate sequences
$usearch -fastx_uniques $merged_fastq_dir/merged_reads_filtered_hmm.fa -minuniquesize 2 -fastaout $merged_fastq_dir/derep_seqs.fa 
# 674248 unique sequences

# Second, sort unique sequences by length
$usearch -sortbylength $merged_fastq_dir/derep_seqs.fa -fastout $merged_fastq_dir/derep_seqs_sorted.fa

# Third, cluster sequences into OTUs at 97% similarity
$usearch -cluster_fast $merged_fastq_dir/derep_seqs.fa -id 0.97 -centroids rep_seqs.fa -uc clusters.uc

#Usearch output
cluster_fast derep_seqs.fa
00:01 172Mb ...... Loading seqs.
00:01 139Mb ...... CPU has 28 cores, defaulting to 10 threads
00:01 172Mb ...... Loading seqs.
00:01 139Mb ...... CPU has 28 cores, defaulting to 10 threads
00:01 152Mb ...... Unique seqs.
00:01 157Mb ...... 260927 seqs, 260927 uniques, 260927 singletons (100.0%)
00:01 157Mb ...... Min size 1, median 1, max 1, avg 1.00
00:01 181Mb 100.0% Convert uniques
00:01 152Mb ...... Unique seqs.
00:01 157Mb ...... 260927 seqs, 260927 uniques, 260927 singletons (100.0%)
00:01 157Mb ...... Min size 1, median 1, max 1, avg 1.00
00:01 181Mb 100.0% Convert uniques
01:02 205Mb 100.0% Unique seqs.
01:02 205Mb 100.0% Writing rep_seqs.fa
01:02 205Mb ......
01:02 205Mb ......       Seqs  260927 (260.9k)
01:02 205Mb ......   Clusters  7489
01:02 205Mb ......   Max size  49002 (49.0k)
01:02 205Mb ......   Avg size  34.8
01:02 205Mb ......   Min size  1
01:02 205Mb ...... Singletons  3067, 1.2% of seqs, 41.0% of clusters
01:02 205Mb ......    Max mem  205Mb
01:02 205Mb ......       Time  01:00
01:02 205Mb ...... Throughput  4348.8 seqs/sec.
01:02 205Mb ......
01:02 205Mb ...... Finished
