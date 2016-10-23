#!/bin/bash
#cd /home/jiangs/jiangs/for_christine
cd /home/christine/lincRNA/testing
# gunzip *

Rscript ./Rscripts/connect_transcripts.R
mkdir ./bed_files
Rscript ./Rscripts/change_and_match.R

mkdir ./neighbors

# need to still sort both lncRNA and coding bed files before using bedtools
sort -k1,1 -k2,2n ./bed_files/lncRNA_matching_transcripts.bed > ./bed_files/lncRNA_matching_transcripts_sorted.bed
sort -k1,1 -k2,2n ./bed_files/coding_matching_transcripts.bed > ./bed_files/coding_matching_transcripts_sorted.bed
# for finding closest neighbors with bedtools
bedtools closest -D ref -a ./bed_files/lncRNA_matching_transcripts_sorted.bed -b ./bed_files/coding_matching_transcripts_sorted.bed > ./neighbors/lncRNA_coding_neighbors
bedtools closest -D ref -a ./bed_files/coding_matching_transcripts_sorted.bed -b ./bed_files/coding_matching_transcripts_sorted.bed > ./neighbors/coding_coding_neighbors

mkdir ./Pearson_graphs
mkdir ./Spearman_graphs
mkdir ./Pearson_values
mkdir ./Spearman_values

Rscript ./Rscripts/calculate_corr_coeff_normalized.R

