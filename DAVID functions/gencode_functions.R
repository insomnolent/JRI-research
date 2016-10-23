# this code is to write gene and transcripts lists into text files for DAVID

# for calculating Pearson Correlation coefficients
# read in both coding_expr and expr_table
codingExpr = read.table("PCG_FPKM", as.is=T, header=T,sep="\t")
noncodingExpr = read.table("lncRNA_FPKM", as.is=T, header=T,sep="\t")
dim(codingExpr) 
#############################
##### for lncRNA:coding #####
#############################
# read in neighbors
neighbors = read.table("./neighbors/lncRNA_coding_neighbors", as.is=T, header=F)
neighbors <- n_lncRNA
dim(neighbors)
# filter out neighbors that are more than 10kb away
neighbors <- neighbors[abs(neighbors[,15]) <= 10000,]
# col for longest overlapping parts
max_right <- ifelse(neighbors[,9]>=neighbors[,2] & neighbors[,15]==0, neighbors[,3]-neighbors[,9], 0)
max_left <- ifelse(neighbors[,9]<=neighbors[,2] & neighbors[,15]==0, neighbors[,10]-neighbors[,2], 0)
max <- max_right+max_left
neighbors <- cbind(neighbors[,1:15], max)
# col for longest nonoverlapping parts of transcripts that overlap
max_right <- ifelse(neighbors[,9]>=neighbors[,2] & neighbors[,15]==0, neighbors[,10]-neighbors[,3], 0)
max_left <- ifelse(neighbors[,9]<=neighbors[,2] & neighbors[,15]==0, neighbors[,2]-neighbors[,9], 0)
max_no_overlap <- max_right+max_left
neighbors <- cbind(neighbors[,1:16], max_no_overlap)
# split col 4 and col 11
library(stringr)
col1 <- str_split_fixed(neighbors[,4], ":", 2)
col2 <- str_split_fixed(neighbors[,11], ":", 2)
neighbors <- data.frame(cbind(neighbors[,1:3],col1,neighbors[,5:10],col2,neighbors[,12:17]))
# pick transcripts with longest overlapping lengths, then longest nonoverlapping lengths
dupes <- neighbors[order(neighbors[,5], -xtfrm(neighbors[,18]),-xtfrm(neighbors[,19])),]
dupes <- dupes[!duplicated(dupes[,5]),]
dim(dupes)
# subset into only gene_IDs and transcript_IDs of neighbors
neighbors <- cbind(dupes[,4:5],dupes[12:13])
# remove rows where transcripts are diff but gene_IDs are the same
neighbors <- neighbors[as.character(neighbors[,1])!=as.character(neighbors[,3]),]
lncRNA_neighbors <- neighbors
write.table(lncRNA_neighbors, "lncRNA_coding_genes.txt",sep="\t",row.names=F,col.names=F,quote=F)
lncRNA_neighbors[,3] <- as.character(lncRNA_neighbors[,3])
test <- str_split_fixed(lncRNA_neighbors[,3], "\\.", 2)
lncRNA_neighbors_genes_only <- test[,1]
write.table(lncRNA_neighbors_genes_only, "lncRNA_coding_genes_only.txt",sep="\t",row.names=F,col.names=F,quote=F)
lncRNA_neighbors[,4] <- as.character(lncRNA_neighbors[,4])
test <- str_split_fixed(lncRNA_neighbors[,4], "\\.", 2)
lncRNA_neighbors_transcript_only <- test[,1]
write.table(lncRNA_neighbors_transcript_only, "lncRNA_coding_transcript_only.txt",sep="\t",row.names=F,col.names=F,quote=F)

#############################
##### for coding:coding #####
#############################
# read in neighbors - has overlap 
n = read.table("./neighbors/coding_coding_neighbors", as.is=T, header=F)
n_coding <- read.table("coding_coding_neighbors_revised", as.is=T, header=F)
neighbors <- n_coding[n_coding$V4!=n_coding$V11,]
dim(neighbors)
# to remove duplicate positions
neighbors <- neighbors[neighbors[,1]==neighbors[,8] & neighbors[,4]!=neighbors[,11],]
neighbors <- neighbors[neighbors[,1]==neighbors[,8] & neighbors[,2]!=neighbors[,9] & neighbors[,3]!=neighbors[,10],]
dim(neighbors)
# if a sequence is completely overlapped by another neighboring coding gene
neighbors <- neighbors[abs(neighbors[,15]) <= 10000,]
dim(neighbors)

# col for longest overlapping parts
max_right <- ifelse(neighbors[,9]>=neighbors[,2] & neighbors[,15]==0, neighbors[,3]-neighbors[,9], 0)
max_left <- ifelse(neighbors[,9]<=neighbors[,2] & neighbors[,15]==0, neighbors[,10]-neighbors[,2], 0)
max <- max_right+max_left
neighbors <- cbind(neighbors[,1:15], max)
# col for longest nonoverlapping parts of transcripts that overlap
max_right <- ifelse(neighbors[,9]>=neighbors[,2] & neighbors[,15]==0, neighbors[,10]-neighbors[,3], 0)
max_left <- ifelse(neighbors[,9]<=neighbors[,2] & neighbors[,15]==0, neighbors[,2]-neighbors[,9], 0)
max_no_overlap <- max_right+max_left
neighbors <- cbind(neighbors[,1:16], max_no_overlap)
# split col 4 and col 11
library(stringr)
col1 <- str_split_fixed(neighbors[,4], ":", 2)
col2 <- str_split_fixed(neighbors[,11], ":", 2)
neighbors <- data.frame(cbind(neighbors[,1:3],col1,neighbors[,5:10],col2,neighbors[,12:17]))
# pick transcripts with longest overlapping lengths, then longest nonoverlapping lengths
dupes <- neighbors[order(neighbors[,5], -xtfrm(neighbors[,18]),-xtfrm(neighbors[,19])),]
dupes <- dupes[!duplicated(dupes[,5]),]
dim(dupes)
# subset into only gene_IDs and transcript_IDs of neighbors
neighbors <- cbind(dupes[,4:5],dupes[12:13])
# remove rows where transcripts are diff but gene_IDs are the same
neighbors <- neighbors[as.character(neighbors[,1])!=as.character(neighbors[,3]),]
# get rid of duplicate neighbors (where there's something like gene1:gene2 and also gene2:gene1)
a <- neighbors[neighbors[,1] %in% neighbors[,3] & neighbors[,3] %in% neighbors[,1],]
a <- a[!duplicated(a[,3]),]
b <- neighbors[!(neighbors[,1] %in% neighbors[,3] & neighbors[,3] %in% neighbors[,1]),]
c <- rbind(b,a)
neighbors <- c
coding_neighbors <- neighbors
write.table(coding_neighbors, "coding_coding_genes.txt",sep="\t",row.names=F,col.names=F,quote=F)

coding_neighbors[,3] <- as.character(coding_neighbors[,3])
test <- str_split_fixed(coding_neighbors[,3], "\\.", 2)
coding_neighbors_genes_only <- test[,1]
write.table(coding_neighbors_genes_only, "coding_coding_genes_only.txt",sep="\t",row.names=F,col.names=F,quote=F)
coding_neighbors[,4] <- as.character(coding_neighbors[,4])
test <- str_split_fixed(coding_neighbors[,4], "\\.", 2)
coding_neighbors_transcript_only <- test[,1]
write.table(coding_neighbors_transcript_only, "coding_coding_transcript_only.txt",sep="\t",row.names=F,col.names=F,quote=F)


# for PCG_PCG_lncRNA

#lncRNA_PCG_PCG <- cbind(lncRNA_left[,1:2], PCG2_right[,1:4])
#PCG_PCG_lncRNA <- cbind(PCG2_left[,3:4], PCG2_left[,1:2],lncRNA_right[,1:2])
write.table(lncRNA_PCG[,3:4], "middle_PCG.txt",sep="\t",row.names=F,col.names=F,quote=F)
# gene_IDs of middle PCGs
c <- as.character(lncRNA_PCG[,3])
test <- str_split_fixed(c, "\\.", 2)
middle_PCG_gene_only <- test[,1]
write.table(middle_PCG_gene_only, "middle_PCG_gene_only.txt",sep="\t",row.names=F,col.names=F,quote=F)
# transcript_IDs of middle PCGs
c <- as.character(lncRNA_PCG[,4])
test <- str_split_fixed(c, "\\.", 2)
middle_PCG_transcript_only <- test[,1]
write.table(middle_PCG_transcript_only, "middle_PCG_transcript_only.txt",sep="\t",row.names=F,col.names=F,quote=F)

