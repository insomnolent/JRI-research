# change it into bed format
gtf <- read.table("gencode_connected_transcripts_nodupes.txt",as.is=T,header=F,sep="\t")
# make new names with both gene_id and transcript_id
a <- gtf[,1]
b <- gtf[,c(4:5)]
c <- paste(gtf[,10],":",gtf[,12],sep="")
d <- gtf[,6:8]
bed <- cbind(a,b,c,d)
write.table(bed, "./bed_files/gencode_connected_transcripts.bed",sep="\t",quote=F,row.names=F,col.names=F)

# all traits (in separate excel sheet)
lncRNA <- read.table("lncRNA_FPKM",as.is=T,header=T,sep="\t")

# find matching locations for each lncRNA
lncRNA_match <- gtf[match(lncRNA[,1], gtf[,12]),]
a <- lncRNA_match[,1]
b <- lncRNA_match[,c(4:5)]
c <- paste(lncRNA_match[,10],":",lncRNA_match[,12],sep="")
d <- lncRNA_match[,6:8]
lncRNA_bed <- cbind(a,b,c,d)
lncRNA_bed <- lncRNA_bed[order(lncRNA_bed[,1]),]
lncRNA_bed <- lncRNA_bed[!is.na(lncRNA_bed[,1]),]
write.table(lncRNA_bed,"./bed_files/lncRNA_matching_transcripts.bed",sep="\t",quote=F,row.names=F,col.names=F)

# all coding genes 
coding <- read.table("PCG_FPKM",as.is=T,header=T,sep="\t")
gtf <- bed

# make sure the gene_id and transcript_id names are in separate columns for the gtf file
library(stringr)
gene_transcript <- str_split_fixed(gtf[,4], ":", 2)
a <- gtf[,1:3]
b <- gene_transcript
c <- gtf[,5:7]
gtf <- cbind(a,b,c)

# find matching locations for each coding
coding_match <- gtf[match(coding[,1], gtf[,5]),]
a <- coding_match[,1:3]
b <- paste(coding_match[,4],":",coding_match[,5],sep="")
c <- coding_match[,6:8]
coding_bed <- cbind(a,b,c)
coding_bed <- coding_bed[order(coding_bed[,1]),]
coding_bed <- coding_bed[!is.na(coding_bed[,1]),]
write.table(coding_bed,"./bed_files/coding_matching_transcripts.bed",sep="\t",quote=F,row.names=F,col.names=F)
