##########      OTU Clustering           ##########
###################################################
#!/bin/bash
# First, dereplicate the sequences
module purge
module load MiniForge3
conda activate

nifH_hmm=/mnt/home/kristybr/NifMAP/Resources/hmm_nuc_1160_nifH.hmm
merged_fastq_dir=/mnt/research/EvansLab/Brandon/Ngrad_2024_nifH/merged_fastq
usearch=/mnt/research/EvansLab/Software/usearch_linux_x86_12.0-beta

# Dereplicate sequences
$usearch -fastx_uniques $merged_fastq_dir/merged_reads_filtered_hmm.fa -minuniquesize 2 -fastaout $merged_fastq_dir/derep_seqs.fa 
# 674248 unique sequences

# Second, sort unique sequences by length
$usearch -sortbylength $merged_fastq_dir/derep_seqs.fa -fastout $merged_fastq_dir/derep_seqs_sorted.fa

# Third, cluster sequences into OTUs at 97% similarity
$usearch -cluster_fast $merged_fastq_dir/derep_seqs.fa -id 0.97 -centroids rep_seqs.fa -uc clusters.uc

