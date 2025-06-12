#!/bin/bash
module purge
module load Miniforge3
conda activate hmmer

nifH_hmm=/mnt/research/EvansLab/Isabella/nifH_tools/nifH_hmms/hmm_nuc_1160_nifH.hmm
merged_fastq_dir=/mnt/research/EvansLab/Isabella/nifH_tools/pellicano_analysis/syncom/syn_merged_fastq/


# Screen merged, filtered reads using HMMER
hmmsearch --domtblout hmmOut1.out $nifH_hmm $merged_fastq_dir/merged_reads_filtered.fa

# Identify acceptable and shit hits 
awk '{print $1}' hmmOut1.out | grep -v "#" > acceptable_hits
# blank reads passed the homolog removal (blank reads removed)
grep ">" $merged_fastq_dir/merged_reads_filtered.fa | grep -v -F -f acceptable_hits > non_nifH_hits
totalUnique=grep ">" $merged_fastq_dir/merged_reads_filtered.fa | wc - l

totalAccepted = `cat acceptable_hits | wc -l`
totalRemoved = `cat non_nifH_hits | wc -l`

# Filter out unacceptable reads from merged_reads_filtered.fa
awk 'BEGIN{FS="\n";RS=">"};NR>1{print(">"$1);for(i=2;i<=NF;i++){printf($i)};print("")}' $merged_fastq_dir/merged_reads_filtered.fa | grep -A 1 -F -f acceptable_hits | grep -v "^\-\-$" > $merged_fastq_dir/merged_reads_filtered_hmm.fa
