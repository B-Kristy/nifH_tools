#!/bin/bash
module purge
module load Miniforge3
conda activate hmmer

nifH_hmm=/mnt/research/EvansLab/Isabella/nifH_tools/nifH_hmms/nifhmm_nuc_1160_nifH.hmm
merged_fastq_dir=/mnt/research/EvansLab/Isabella/merged_fastq


# Screen merged, filtered reads using HMMER
hmmsearch --domtblout hmmOut1.out $nifH_hmm $merged_fastq_dir/merged_reads_filtered.fa