gtf <- read.table("gencode_refseq_lncRNAdb_V5.gtf",as.is=T,header=F,sep="\t")
# split last description col into mult cols
library(stringr)
col <- str_split_fixed(gtf[,9], "; ", 3)
gtf <- cbind(gtf[,1:8],col[,1:2])
# filter out everything except exons
exon <- gtf[gtf[,3] == "exon",]
# filter out things that don't have transcript IDs
a <- exon[!grepl("transcript_id",exon[,10]),]
# split gene_id and transcript_id into their own column
gene_id <- str_split_fixed(exon[,9], " ", 2)
transcript_id <- str_split_fixed(exon[,10], " ", 2)
exon <- cbind(exon[,1:8],gene_id,transcript_id)
# sort exons by chromosomes
transcript <- exon[order(exon[,1]),]

#########################
## connect transcripts ##
#########################
# separate based on database source: HAVANA, ENSEMBL, RefSeq, lncRNAdb
sep <- unique(transcript[,2])
havana <- transcript[grep("HAVANA", transcript[,2]),]
ensembl <- transcript[grep("ENSEMBL", transcript[,2]),]
refseq <- transcript[grep("RefSeq", transcript[,2]),]
lncRNAdb <- transcript[grep("lncRNAdb", transcript[,2]),]
# sort each subset by position
havana <- havana[order(havana[,4]),]
ensembl <- ensembl[order(ensembl[,4]),]
refseq <- refseq[order(refseq[,4]),]
lncRNAdb <- lncRNAdb[order(lncRNAdb[,4]),]
# filter out start positions of first unique transcript id
havana_a <- havana[!duplicated(havana[,12]),]
havana_start <- havana_a[,4]
ensembl_a <- ensembl[!duplicated(ensembl[,12]),]
ensembl_start <- ensembl_a[,4]
refseq_a <- refseq[!duplicated(refseq[,12]),]
refseq_start <- refseq_a[,4]
lncRNAdb_a <- lncRNAdb[!duplicated(lncRNAdb[,12]),]
lncRNAdb_start <- lncRNAdb_a[,4]
# filter out end positions of first unique transcript id
b <- havana[!rev(duplicated(rev(havana[,12]))),]
b <- havana[!duplicated(havana[,12], fromLast=T), ]
havana_end <- b[,5]
b <- ensembl[!rev(duplicated(rev(ensembl[,12]))),]
b <- ensembl[!duplicated(ensembl[,12], fromLast=T), ]
ensembl_end <- b[,5]
b <- refseq[!rev(duplicated(rev(refseq[,12]))),]
b <- refseq[!duplicated(refseq[,12], fromLast=T), ]
refseq_end <- b[,5]
b <- lncRNAdb[!rev(duplicated(rev(lncRNAdb[,12]))),]
b <- lncRNAdb[!duplicated(lncRNAdb[,12], fromLast=T), ]
lncRNAdb_end <- b[,5]
# connect exons into transcripts
havana_connected <- cbind(havana_a[,c(1:3)],havana_start, havana_end, havana_a[,c(6:12)])
ensembl_connected <- cbind(ensembl_a[,c(1:3)],ensembl_start, ensembl_end, ensembl_a[,c(6:12)])
refseq_connected <- cbind(refseq_a[,c(1:3)],refseq_start, refseq_end, refseq_a[,c(6:12)])
lncRNAdb_connected <- cbind(lncRNAdb_a[,c(1:3)],lncRNAdb_start, lncRNAdb_end, lncRNAdb_a[,c(6:12)])
# sort each connected subset by chromosomes
havana_connected <- havana_connected[order(havana_connected[,1]),]
ensembl_connected <- ensembl_connected[order(ensembl_connected[,1]),]
refseq_connected <- refseq_connected[order(refseq_connected[,1]),]
lncRNAdb_connected <- lncRNAdb_connected[order(lncRNAdb_connected[,1]),]
# need to check that all values in end are greater than start - with sum(x > 0)
# combine all the different transcripts
connected <- rbind(as.matrix(havana_connected),as.matrix(ensembl_connected),as.matrix(refseq_connected),as.matrix(lncRNAdb_connected))
connected[,12] <- sub(";","",connected[,12])
write.table(connected, "gencode_connected_transcripts_nodupes.txt",sep="\t",quote=F,row.names=F,col.names=F)
