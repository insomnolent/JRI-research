# this is for checking the overlap between gene lists and transcript lists for coding neighbors

# overlap between lincRNA neighbors gene and middle gene
length(intersect(lncRNA_neighbors_genes_only,middle_PCG_gene_only))
# overlap between lincRNA neighbors transcript and middle transcript
length(intersect(lncRNA_neighbors_transcript_only,middle_PCG_transcript_only))
# check for overlap later