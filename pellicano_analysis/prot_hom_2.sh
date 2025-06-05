#!/bin/bash


cd /mnt/research/EvansLab/Isabella/merged_fastq/
cat nifH_corr_prot.fasta | sed 's/ //g' |  awk 'BEGIN{RS=">";FS="\n"};NR>1{printf $1"_\t";for(i=2;i<=NF;i++){printf($i)}print("")}'
cd /mnt/research/EvansLab/Isabella/nifH_tools/pellicano_analysis
grep -F -f protein_homologs | awk '{gsub("_$","",$1);print(">"$1"\n"$2)}' > nifH_rej_prot.fasta