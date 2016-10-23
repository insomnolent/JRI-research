#!/bin/bash

# change .bed to .bedg
cd /home/christine/
bedtools map -a ./lincRNA/testing/bed_files/lncRNA_matching_transcripts_sorted.bed -b ./phylop/out.bedg -c 4 -o mean > ./phylop/lncRNA_scores

bedtools map -a ./lincRNA/testing/bed_files/coding_matching_transcripts_sorted.bed -b ./phylop/out.bedg -c 4 -o mean > ./phylop/coding_scores