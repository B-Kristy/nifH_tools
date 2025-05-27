#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########

#SBATCH --time=40:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=32                   # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=32                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=16          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=128G          # memory required per allocated CPU (or core)
#SBATCH --job-name merge_seq       # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kristybr@msu.edu

########## Command Lines for Job Running ##########

##########      Rename fastq files       ##########
###################################################


for f in *.fastq.gz; do
  newname=$(echo "$f" | sed -E 's/S[0-9]+_L001_//')
  mv "$f" "$newname"
done

for f in *.fastq.gz; do
  # Extract the first R[0-9]+ match
  first_r=$(echo "$f" | grep -oE 'R[0-9]+' | head -n1)

  # Count total R[0-9]+ matches
  rcount=$(echo "$f" | grep -oE 'R[0-9]+' | wc -l)

  if [[ $rcount -ge 2 ]]; then
    # Build replacement string
    repx="Rep${first_r:1}"  # strip the 'R' and keep the number
    # Replace only the first R[0-9]+ with RepxX
    newname=$(echo "$f" | sed -E "s/${first_r}/${repx}/")
    mv "$f" "$newname"
  fi
done


########## Assemble Paired-End Sequences ##########
###################################################
#!/bin/bash
fastq_dir=/mnt/research/EvansLab/Brandon/Ngrad_2024_nifH/raw_fastq
merged_fastq_dir=/mnt/research/EvansLab/Brandon/Ngrad_2024_nifH/merged_fastq
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
cat *merged_reads* > merged_reads.fastq
mv merged_reads.fastq $merged_fastq_dir

# Remove intermediate merge files
rm -rf $merged_fastq_dir/temp
mkdir $merged_fastq_dir/temp

# Count the number of total reads before filtering 
wc -l merged_reads.fastq | awk '{print $1/4}' # 5,060,842 reads total (~22,000 reads per sample)

##########    Filter assembled reads     ##########
###################################################
#!/bin/bash
merged_fastq_dir=/mnt/research/EvansLab/Brandon/Ngrad_2024_nifH/merged_fastq
usearch=/mnt/research/EvansLab/Software/usearch_linux_x86_12.0-beta

cd $merged_fastq_dir
$usearch -fastq_filter merged_reads.fastq -fastq_trunclen 380 -fastq_maxee 1.0 -fastaout merged_reads_filtered.fa

# Filtered reads result
02:34 2.2Mb ......    5060842  Reads (5.1M)
02:34 2.2Mb ......      15356  Discarded reads length < 380
02:34 2.2Mb ......     331303  Discarded reads with expected errs > 1.00
02:34 2.2Mb ......    4714183  Filtered reads (4.7M, 93.2%)
02:34 2.2Mb ...... Finished

# Overall, 4,714,183 reads were retains (93.2% of total dataset)

##########       Homolog Removal         ##########
###################################################
#!/bin/bash
module purge
module load Conda/3

nifH_hmm=/mnt/home/kristybr/NifMAP/Resources/hmm_nuc_1160_nifH.hmm
merged_fastq_dir=/mnt/research/EvansLab/Brandon/Ngrad_2024_nifH/merged_fastq


# Screen merged, filtered reads using HMMER
hmmsearch --domtblout hmmOut1.out $nifH_hmm $merged_fastq_dir/merged_reads_filtered.fa

# Identify acceptable and shit hits 
awk '{print $1}' hmmOut1.out | grep -v "#" > acceptable_hits
# 4,712,853 reads passed the homolog removal (1,332 reads removed)
grep ">" $merged_fastq_dir/merged_reads_filtered.fa | grep -v -F -f acceptable_hits > non_nifH_hits
totalUnique=grep ">" $merged_fastq_dir/merged_reads_filtered.fa | wc - l

totalAccepted = `cat acceptable_hits | wc -l`
totalRemoved = `cat non_nifH_hits | wc -l`

# Filter out unacceptable reads from merged_reads_filtered.fa
awk 'BEGIN{FS="\n";RS=">"};NR>1{print(">"$1);for(i=2;i<=NF;i++){printf($i)};print("")}' $merged_fastq_dir/merged_reads_filtered.fa | grep -A 1 -F -f acceptable_hits | grep -v "^\-\-$" > $merged_fastq_dir/merged_reads_filtered_hmm.fa

##########      OTU Clustering           ##########
###################################################
#!/bin/bash
# First, dereplicate the sequences
module purge
module load Conda/3

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

# 15,997 OTUs were generated at 97% similarity 

# Usearch output
cluster_fast derep_seqs.fa
00:03 390Mb ...... Loading seqs.
00:03 356Mb ...... CPU has 128 cores, defaulting to 10 threads
00:03 390Mb ...... Unique seqs.
00:03 404Mb ...... 674248 seqs, 674248 uniques, 674248 singletons (100.0%)
00:03 404Mb ...... Min size 1, median 1, max 1, avg 1.00
00:04 473Mb 100.0% Convert uniques
04:00 517Mb 100.0% Unique seqs.
04:00 517Mb 100.0% Writing rep_seqs.fa
04:00 517Mb ...... 
04:00 517Mb ......       Seqs  674248 (674.2k)
04:00 517Mb ......   Clusters  15997 (16.0k)
04:00 517Mb ......   Max size  94838 (94.8k)
04:00 517Mb ......   Avg size  42.1
04:00 517Mb ......   Min size  1
04:00 517Mb ...... Singletons  6957, 1.0% of seqs, 43.5% of clusters
04:00 517Mb ......    Max mem  517Mb
04:00 517Mb ......       Time  03:56
04:00 517Mb ...... Throughput  2857.0 seqs/sec.
04:00 517Mb ...... 
04:00 517Mb ...... Finished

########## Protein Sequence Conversion   ##########
###################################################
#!/bin/bash
merged_fastq_dir=/mnt/research/EvansLab/Brandon/Ngrad_2024_nifH/merged_fastq
framebot=/mnt/home/kristybr/anaconda3/pkgs/rdptools-2.0.2-1/bin/FrameBot
nifH_prot_ref=/mnt/home/kristybr/Framebot-master/refset/nifh_prot_ref.fasta

# Translate OTU representative sequences into protein sequences
cd $merged_fastq_dir
$framebot framebot -N -l 30 -i 0.4 -o $merged_fastq_dir $nifH_prot_ref $merged_fastq_dir/rep_seqs.fa


########## Protein Homolog Removal       ##########
###################################################
#!/bin/bash

module purge
module load Conda/3

merged_fastq_dir=/mnt/research/EvansLab/Brandon/Ngrad_2024_nifH/merged_fastq
prot_hmm_dir=/mnt/home/kristybr/NifMAP/Resources/
# Corrected AA sequences are in _corr_prot.fasta
# Screen with hmm to identify all hits

# First - press the homolog hmm prior to running hmmscan - code will not run unless you do this
hmmpress $prot_hmm_dir/nifH_chlL_bchX.hmm

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
totalOTUs=grep ">" $merged_fastq_dir/merged_fastq_corr_nucl.fasta | wc -l # 15,992 OTUs
totalAccepted=cat acceptable_hits | wc -l # 6,030 OTUs
totalRemoved=cat protein_homologs | wc -l # 9,962 OTUs removed as partial nifH or homologous protein sequences 


#mv ${RESULTSFOLDER}/${inputReads} ${WORKFOLDER}/${inputReads}
cat $merged_fastq_dir/merged_fastq_corr_nucl.fasta | sed 's/ //g' | awk 'BEGIN{RS=">";FS="\n"};NR>1{printf $1"_\t";for(i=2;i<=NF;i++){printf($i)}print("")}' | grep -F -f acceptable_hits | awk '{gsub("_$","",$1);print(">"$1"\n"$2)}' > nifH_corr_nucl_only_nifH.fasta
cat $merged_fastq_dir/merged_fastq_corr_prot.fasta | sed 's/ //g' |  awk 'BEGIN{RS=">";FS="\n"};NR>1{printf $1"_\t";for(i=2;i<=NF;i++){printf($i)}print("")}' | grep -F -f acceptable_hits | awk '{gsub("_$","",$1);print(">"$1"\n"$2)}' > nifH_corr_prot.fasta
cat nifH_corr_prot.fasta | sed 's/ //g' |  awk 'BEGIN{RS=">";FS="\n"};NR>1{printf $1"_\t";for(i=2;i<=NF;i++){printf($i)}print("")}' | grep -F -f protein_homologs | awk '{gsub("_$","",$1);print(">"$1"\n"$2)}' > nifH_rej_prot.fasta

##########    OTU Table Construction     ##########
###################################################
#!/bin/bash

module purge
module load Conda/3

merged_fastq_dir=/mnt/research/EvansLab/Brandon/Ngrad_2024_nifH/merged_fastq
usearch=/mnt/research/EvansLab/Software/usearch_linux_x86_12.0-beta


# Make an OTU table based on filter corrected seqs
$usearch -otutab $merged_fastq_dir/merged_reads_filtered_hmm.fa -otus $merged_fastq_dir/nifH_corr_nucl_only_nifH.fasta -otutabout otu_table.txt -mapout map.txt -threads 32

# OTU Table Generation Output 
00:00 7.2Mb ...... Loading seqs.
00:00 4.5Mb 100.0% Masking (fastnucleo)
00:00 5.5Mb 100.0% Udb seqs
00:00 5.5Mb 100.0% Alloc rows
00:00  14Mb 100.0% Build udb index
00:00  14Mb ...... Starting 32 threads
03:04 111Mb 100.0% Mapping 6030 OTUs, 228 samples, 1.8M assigned (37.4%)  
03:04 111Mb ...... 1760408 / 4712851 mapped to OTUs (37.4%)
03:04 111Mb ...... Finished

##########   nifH Feature Classifier     ##########
###################################################
#!/bin/bash

module purge
module load Conda/3
conda activate qiime2-amplicon-2024.10 

nifH_dada2_dir=/mnt/home/kristybr/nifHdada2-master
merged_fastq_dir=/mnt/research/EvansLab/Brandon/Ngrad_2024_nifH/merged_fastq

cd $merged_fastq_dir
# Make taxonomy file
touch nifH_dada2_v2.0.5_tax.tsv
echo -e "Feature ID\tTaxon" >> nifH_dada2_v2.0.5_tax.tsv
tail -n +2 $nifH_dada2_dir/nifH_dada2_phylum_v2.0.5.csv | cut -f 1,5 -d ',' | sed 's/,/\t/' >> nifH_dada2_v2.0.5_tax.tsv

# Make fasta file 
tail -n +2 $nifH_dada2_dir/nifH_dada2_phylum_v2.0.5.csv | cut -f 1,2 -d ',' | sed 's/^/>/' | sed 's/,/\n/'  > nifH_dada2_v2.0.5_seqs.fasta

# Import into QIIME2
qiime tools import --input-path nifH_dada2_v2.0.5_tax.tsv --type 'FeatureData[Taxonomy]' --output-path nifH_dada2_v2.0.5_tax.qza
qiime tools import --input-path nifH_dada2_v2.0.5_seqs.fasta --type 'FeatureData[Sequence]' --output-path nifH_dada2_v2.0.5_seqs.qza

# Train feature classifier
qiime feature-classifier fit-classifier-naive-bayes  --i-reference-reads nifH_dada2_v2.0.5_seqs.qza --i-reference-taxonomy nifH_dada2_v2.0.5_tax.qza --o-classifier nifH_dada2_v2.0.5_classifier

##########Import filtered sequences into QIIME####
###################################################
#!/bin/bash

merged_fastq_dir=/mnt/research/EvansLab/Brandon/Ngrad_2024_nifH/merged_fastq
cd $merged_fastq_dir

qiime tools import --input-path nifH_corr_nucl_only_nifH.fasta --output-path nifH_sequences.qza --type 'FeatureData[Sequence]'

##########      Classify Sequences       ##########
###################################################
#!/bin/bash

module purge
module load Conda/3
conda activate qiime2-amplicon-2024.10 

nifH_dada2_dir=/mnt/home/kristybr/nifHdada2-master
merged_fastq_dir=/mnt/research/EvansLab/Brandon/Ngrad_2024_nifH/merged_fastq

cd $merged_fastq_dir

qiime feature-classifier classify-sklearn --i-classifier nifH_dada2_v2.0.5_classifier.qza --i-reads nifH_sequences.qza --o-classification nifH_taxonomy.qza

qiime metadata tabulate --m-input-file nifH_taxonomy.qza  --o-visualization nifH_taxonomy.qzv


##########     Re-Do Taxonomy from NFixDB       ###
###################################################
#!/bin/bash

# Lisa Tiemann compiled the NFixDB and trained the feature classifier  

module purge
module load Conda/3
conda activate qiime2-amplicon-2024.10 

nfix_db_dir=/mnt/home/kristybr/NFixDB
merged_fastq_dir=/mnt/research/EvansLab/Brandon/Ngrad_2024_nifH/merged_fastq

cd $merged_fastq_dir

qiime feature-classifier classify-sklearn --i-classifier $nfix_db_dir/nifH_classifier.qza --i-reads nifH_sequences.qza --o-classification nifH_taxonomy.qza
qiime metadata tabulate --m-input-file nifH_taxonomy.qza  --o-visualization nifH_taxonomy.qzv
