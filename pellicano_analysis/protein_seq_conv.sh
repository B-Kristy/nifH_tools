#!/bin/bash
module purge
module load Miniforge3
conda activate framebot

merged_fastq_dir=/mnt/research/EvansLab/Isabella/merged_fastq
nifH_prot_ref=/mnt/research/EvansLab/Isabella/nifH_tools/nifh_prot_ref.fasta
framebot=/mnt/home/f0111521/RDPTools/

cd $framebot
ant -f Framebot/build.xml jar


# Translate OTU representative sequences into protein sequences
cd $merged_fastq_dir
$framebot FrameBot framebot -N -l 30 -i 0.4 -o $merged_fastq_dir $nifH_prot_ref $merged_fastq_dir/rep_seqs.fa
