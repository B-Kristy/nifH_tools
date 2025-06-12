#!/bin/bash
fastq_dir=/mnt/research/EvansLab/Isabella/nifH_tools/pellicano_analysis/mock_c07
merged_fastq_dir=/mnt/research/EvansLab/Isabella/merged_fastq/
usearch=/mnt/research/EvansLab/Software/usearch_linux_x86_12.0-beta

cd $fastq_dir
# Merge paired-end fastq sequences using usearch
for i in *R1_001*.fastq; do
  $usearch
  -fastq_mergepairs $i \
  -fastqout $merged_fastq_dir/merged_reads_${i}\
  -sample $i
done

# Move intermediate files into temp storage 
mkdir $merged_fastq_dir/temp
mv *.fastq ./temp

# Concatenate all merged reads into one file
cd $merged_fastq_dir/temp
cat *merged_reads* > c07_merged_reads.fastq
mv c07_merged_reads.fastq $merged_fastq_dir

# Remove intermediate merge files
rm -rf $merged_fastq_dir/temp
mkdir $merged_fastq_dir/temp

# Count the number of total reads before filtering 
wc -l c07_merged_reads.fastq | awk '{print $1/4}' # 5,060,842 reads total (~22,000 reads per sample)