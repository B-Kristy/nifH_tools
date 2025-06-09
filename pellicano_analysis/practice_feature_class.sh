#!/bin/bash

module purge
module load Miniforge3
conda activate qiime2

#get reference data set and import as qiime artifacts
cd /mnt/research/EvansLab/Isabella/nifH_tools/pellicano_analysis

wget -O "85_otu_taxonomy.txt" "https://data.qiime2.org/2024.10/tutorials/training-feature-classifiers/85_otu_taxonomy.txt"

qiime tools import --type 'FeatureData[Sequence]' --input-path 85_otus.fasta --output-path 85_otus.qza
qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path 85_otu_taxonomy.txt --output-path ref-taxonomy.qza

#extract reference reads##
qiime feature-classifier extract-reads --i-sequences 85_otus.qza --p-f-primer GTGCCAGCMGCCGCGGTAA --p-r-primer GGACTACHVGGGTWTCTAAT --p-trunc-len 120 --p-min-length 100 --p-max-length 400 --o-reads ref-seqs.qza

#train the classifier##
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ref-seqs.qza --i-reference-taxonomy ref-taxonomy.qza --o-classifier classifier.qza

#test the classifier##
wget -O "rep-seqs.qza" "https://data.qiime2.org/2024.10/tutorials/training-feature-classifiers/rep-seqs.qza"

qiime feature-classifier classify-sklearn --i-classifier classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza

qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv