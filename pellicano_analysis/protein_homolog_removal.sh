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

grep "nifH" assignments.txt | awk '{print $2}' | sort > acceptable_hits
grep ">" $merged_fastq_dir/merged_fastq_corr_nucl.fasta | awk '{print $1}' | grep -v -F -f acceptable_hits | sed 's/>//'> protein_homologs
totalOTUs=grep ">" $merged_fastq_dir/merged_fastq_corr_nucl.fasta | wc -l 
totalAccepted=cat acceptable_hits | wc -l 
totalRemoved=cat protein_homologs | wc -l 


#mv ${RESULTSFOLDER}/${inputReads} ${WORKFOLDER}/${inputReads}
cat $merged_fastq_dir/merged_fastq_corr_nucl.fasta | sed 's/ //g' | awk 'BEGIN{RS=">";FS="\n"};NR>1{printf $1"_\t";for(i=2;i<=NF;i++){printf($i)}print("")}' | grep -F -f acceptable_hits | awk '{gsub("_$","",$1);print(">"$1"\n"$2)}' > nifH_corr_nucl_only_nifH.fasta
cat $merged_fastq_dir/merged_fastq_corr_prot.fasta | sed 's/ //g' |  awk 'BEGIN{RS=">";FS="\n"};NR>1{printf $1"_\t";for(i=2;i<=NF;i++){printf($i)}print("")}' | grep -F -f acceptable_hits | awk '{gsub("_$","",$1);print(">"$1"\n"$2)}' > nifH_corr_prot.fasta
cat nifH_corr_prot.fasta | sed 's/ //g' |  awk 'BEGIN{RS=">";FS="\n"};NR>1{printf $1"_\t";for(i=2;i<=NF;i++){printf($i)}print("")}' | grep -F -f protein_homologs | awk '{gsub("_$","",$1);print(">"$1"\n"$2)}' > nifH_rej_prot.fasta