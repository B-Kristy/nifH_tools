#!/bin/bash

module purge
module load Miniforge3
conda activate hmmer

merged_fastq_dir=/mnt/research/EvansLab/Isabella/merged_fastq/
prot_hmm_dir=/mnt/research/EvansLab/Isabella/nifH_tools/nifH_hmms/
# Corrected AA sequences are in _corr_prot.fasta
# Screen with hmm to identify all hits

# hmmscan param
eVal_chL=1e-50
score_chL=150
eVal_bchX=1e-50
score_bchX=150

# Now run hmmscan
hmmscan --cpu 32  --domtblout hmmOut2.out $prot_hmm_dir/nifH_chlL_bchX.hmm $merged_fastq_dir/merged_fastq_corr_prot.fasta &> /dev/null

# Filter out homologous sequences
cat hmmOut2.out | awk 'NR>3{if($8>bitarray[$4]){bitarray[$4]=$8;outArray[$4]=$1"\t"$4}}END{for(entry in outArray){print outArray[entry]}}' > assignments.txt