# for calculating Pearson Correlation coefficients
# read in both coding_expr and expr_table
setwd("/Users/christinesun/Documents/Research/JRI/gencode/coding_coding")
codingExpr = read.table("PCG_FPKM", as.is=T, header=T,sep="\t")
noncodingExpr = read.table("lncRNA_FPKM", as.is=T, header=T,sep="\t")
n_lncRNA = read.table("./lncRNA_coding_neighbors_no_overlap", as.is=T, header=F)
n_coding <- read.table("./coding_coding_neighbors_no_overlap", as.is=T, header=F)
#############################
##### for lncRNA:coding #####
#############################
# read in neighbors
neighbors <- n_lncRNA
dim(neighbors)
neighbors <- neighbors[abs(neighbors[,15]) <= 500000,]
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
neighbors <- dupes[as.character(dupes[,4])!=as.character(dupes[,12]),]
lncRNA_distance = abs(neighbors[,17])
fn <- ecdf(lncRNA_distance)
plot(fn, main="Noncoding neighbor distances", xlab="Distance to closest coding neighbors (KB)", ylab="Cumulative proportion") 

#############################
##### for coding:coding #####
#############################
# read in neighbors - has overlap 
neighbors <- n_coding[n_coding$V4!=n_coding$V11,]
dim(neighbors)
# to remove duplicate positions
neighbors <- neighbors[neighbors[,1]==neighbors[,8] & neighbors[,4]!=neighbors[,11],]
neighbors <- neighbors[neighbors[,1]==neighbors[,8] & neighbors[,2]!=neighbors[,9] & neighbors[,3]!=neighbors[,10],]
dim(neighbors)
# if a sequence is completely overlapped by another neighboring coding gene
neighbors <- neighbors[abs(neighbors[,15]) <= 500000,]
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
neighbors <- dupes
neighbors <- neighbors[as.character(neighbors[,4])!=as.character(neighbors[,12]),]
coding_distance = abs(neighbors[,17])

fn <- ecdf(coding_distance)
plot(fn, main="Coding neighbor distances", xlab="Distance to closest coding neighbors (KB)", ylab="Cumulative proportion") 


##################################
## create plot with all 3 lines ##
##################################
# with Pearson values
lncRNA_distance <- data.frame(lncRNA_distance)
lncRNA_distance$name <- 'lncRNA_neighbors'
coding_distance <- data.frame(coding_distance)
coding_distance$name <- 'coding_neighbors'
colnames(lncRNA_distance) <- c(coding_distance, name) 
neighbors_distance <- rbind("lncRNA_distance","coding_distance")

sample1 <- lncRNA_distance[,1]
sample2 <- coding_distance[,1]

plot(ecdf(sample1), verticals=TRUE, do.p=FALSE, main="Neighbor distances", xlab="Distance to closest coding neighbors (bases)",ylab="Cumulative Percent",lty="solid")
lines(ecdf(sample2), verticals=TRUE, do.p=FALSE,col.h="red", col.v="red",lty="solid")
legend(50000,.4,c("lncRNA","coding"),col=c("black","red"),lty=c("solid","solid"), pt.cex=2)
