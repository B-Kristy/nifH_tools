#!/bin/bash

module purge
module load Miniforge3
conda activate qiime2

#Import Buckley Data 
qiime tools import --type 'FeatureData[Sequence]' --input-path Buckley_nifH_seqs_cleaned_qiime2format.fasta --output-path Buckley_nifH_seqs_cleaned_qiime2format.qza
qiime tools import --type 'FeatureData[Taxonomy]' --input-path Buckley_nifH_taxonomy_cleaned_qiime2format.tsv --output-path Buckley_nifH_taxonomy_cleaned_qiime2format.qza

#Train Buckley classifier
qiime feature-classifier fit-classifier-naive-bayes  --i-reference-reads Buckley_nifH_seqs_cleaned_qiime2format.qza --i-reference-taxonomy Buckley_nifH_taxonomy_cleaned_qiime2format.qza --o-classifier Buckley_nifH_classifier.qza

#Import NFixDB data
qiime tools import --type 'FeatureData[Sequence]' --input-path NfixDB_nifH_seqs_qiime2format.fasta --output-path NfixDB_nifH_seqs_qiime2format.qza
qiime tools import --type 'FeatureData[Taxonomy]' --input-path NFixDB_nifH_taxonomy_qiime2format.tsv --output-path NFixDB_nifH_taxonomy_qiime2format.qza

#Train NFixDB classifier
qiime feature-classifier fit-classifier-naive-bayes  --i-reference-reads NfixDB_nifH_seqs_qiime2format.qza --i-reference-taxonomy NFixDB_nifH_taxonomy_qiime2format.qza --o-classifier NFixDB_nifH_classifier.qza

#Import Zehr Data
qiime tools import --type 'FeatureData[Sequence]' --input-path Zehr_nifH_seqs_qiime2format.fasta --output-path Zehr_nifH_seqs_qiime2format.qza
qiime tools import --type 'FeatureData[Taxonomy]' --input-path Zehr_nifH_taxonomy_qiime2format.tsv --output-path Zehr_nifH_taxonomy_qiime2format.qza

#Train Zehr classifier
qiime feature-classifier fit-classifier-naive-bayes  --i-reference-reads Zehr_nifH_seqs_qiime2format.qza --i-reference-taxonomy Zehr_nifH_taxonomy_qiime2format.qza --o-classifier Zehr_nifH_classifier.qza

#import actual data
qiime tools import --type 'FeatureData[Sequence]' --input-path nifH_corr_nucl_only_nifH_upper_.fasta --output-path nifH_sequences.qza 


###CLASSIFY READS###
bigdir=/mnt/research/EvansLab/Isabella/nifH_tools/nifH_databases
buckdir=/mnt/research/EvansLab/Isabella/nifH_tools/nifH_databases/Gaby_Buckley
nfixdir=/mnt/research/EvansLab/Isabella/nifH_tools/nifH_databases/NFixDB
zehrdir=/mnt/research/EvansLab/Isabella/nifH_tools/nifH_databases/Zehr

#classify using Buckley
qiime feature-classifier classify-sklearn --i-classifier $buckdir/Buckley_nifH_classifier.qza --i-reads $bigdir/nifH_sequences.qza --o-classification nifH_taxonomy_buckley.qza
qiime metadata tabulate --m-input-file nifH_taxonomy_buckley.qza  --o-visualization nifH_taxonomy_buckley.qzv

#classify using NFixDB
qiime feature-classifier classify-sklearn --i-classifier $nfixdir/NFixDB_nifH_classifier.qza --i-reads $bigdir/nifH_sequences.qza --o-classification nifH_taxonomy_nfixdb.qza
qiime metadata tabulate --m-input-file nifH_taxonomy_nfixdb.qza  --o-visualization nifH_taxonomy_nfixdb.qzv

#classify using Zehr
qiime feature-classifier classify-sklearn --i-classifier $zehrdir/Zehr_nifH_classifier.qza --i-reads $bigdir/nifH_sequences.qza --o-classification nifH_taxonomy_zehr.qza
qiime metadata tabulate --m-input-file nifH_taxonomy_zehr.qza  --o-visualization nifH_taxonomy_zehr.qzv
