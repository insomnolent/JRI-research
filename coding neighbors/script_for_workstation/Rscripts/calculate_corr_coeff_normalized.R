# for calculating Pearson Correlation coefficients
# read in both coding_expr and expr_table
codingExpr = read.table("PCG_FPKM", as.is=T, header=T,sep="\t")
dim(codingExpr)
noncodingExpr = read.table("lncRNA_FPKM", as.is=T, header=T,sep="\t")
dim(noncodingExpr)
n_lncRNA = read.table("./neighbors/lncRNA_coding_neighbors", as.is=T, header=F)
n_coding = read.table("./neighbors/coding_coding_neighbors", as.is=T, header=F)
#############################
##### for lncRNA:coding #####
#############################
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

# with duplicate gene_IDS included
coding1 <- noncodingExpr[match(neighbors[,2], noncodingExpr[,1]),]
# to normalize lncRNA values
pname <- coding1[,1]
pexpr <- coding1[,2:length(coding1)]
plog2 <- log2(pexpr + 1)
psum <- sum(plog2)
pnorm <- plog2 / psum
coding1 <- cbind(pname,pnorm)
# for PCG genes
coding2 <- codingExpr[match(neighbors[,4], codingExpr[,1]),]
# to normalize PCG values
pname <- coding2[,1]
pexpr <- coding2[,2:length(coding2)]
plog2 <- log2(pexpr + 1)
psum <- sum(plog2)
pnorm <- plog2 / psum
coding2 <- cbind(pname, pnorm)
# continue with analysis
coding <- cbind(coding1, coding2)
coding <- coding[!is.na(coding[,101]),]
coding <- coding[!is.na(coding[,1]),]
coding <- as.matrix(coding)
values <- vector(mode="numeric", length=length(coding[,1]))
for (i in 1:length(coding[,1])){
  values[i] <- cor(as.numeric(coding[i,2:100]),as.numeric(coding[i,102:200]),method="pearson")
  if (i%%1000==0) print(i)
}
# remove NAs from values
values <- values[!is.na(values)]
write.table(values,"./Pearson_values/lncRNA_coding_normalized.txt",sep="\t",quote=F,row.names=F,col.names=F)

# plot densities of the values 
library(ggplot2)
X <- values
data <- data.frame(X)
pdf("./Pearson_graphs/lncRNA_coding_normalized_pearson.pdf")
ggplot(data, aes(x = X)) + geom_density() + labs(title = "Neighbors", x="Spearman Correlation", y="Density")
dev.off()
pdf("./Pearson_graphs/lncRNA_coding_normalized_histogram_pearson.pdf")
hist(values)
dev.off()

# calc Spearman correlations too
values <- vector(mode="numeric", length=length(coding[,1]))
for (i in 1:length(coding[,1])){
    values[i] <- cor(as.numeric(coding[i,2:100]),as.numeric(coding[i,102:200]),method="spearman")
    if (i%%1000==0) print(i)
}
# remove NAs from values
values <- values[!is.na(values)]
write.table(values,"./Spearman_values/lncRNA_coding_normalized_spearman.txt",sep="\t",quote=F,row.names=F,col.names=F)

pdf("./Spearman_graphs/lncRNA_coding_normalized_spearman.pdf")
ggplot(data, aes(x = X)) + geom_density() + labs(title = "Neighbors", x="Spearman Correlation", y="Density")
dev.off()
pdf("./Spearman_graphs/lncRNA_coding_normalized_histogram_spearman.pdf")
hist(values)
dev.off()

#############################
##### for coding:coding #####
#############################
neighbors <- n_coding
dim(neighbors)
neighbors <- neighbors[neighbors$V4!=neighbors$V11,]
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
# neighbors without duplicate gene_IDs
coding1 <- codingExpr[match(neighbors[,2], codingExpr[,1]),]
# to normalize lncRNA values
pname <- coding1[,1]
pexpr <- coding1[,2:length(coding1)]
plog2 <- log2(pexpr + 1)
psum <- sum(plog2)
pnorm <- plog2 / psum
coding1 <- cbind(pname,pnorm)
# PCG's neighbors
coding2 <- codingExpr[match(neighbors[,4], codingExpr[,1]),]
pname <- coding2[,1]
pexpr <- coding2[,2:length(coding2)]
plog2 <- log2(pexpr + 1)
psum <- sum(plog2)
pnorm <- plog2 / psum
coding2 <- cbind(pname, pnorm)
# rest of analysis
coding <- cbind(coding1, coding2)
coding <- coding[!is.na(coding[,101]),]
coding <- coding[!is.na(coding[,1]),]
coding <- as.matrix(coding)
values <- vector(mode="numeric", length=length(coding[,1]))
for (i in 1:length(coding[,1])){
  values[i] <- cor(as.numeric(coding[i,2:100]),as.numeric(coding[i,102:200]),method="pearson")
  if (i%%1000==0) print(i)
}
values <- values[!is.na(values)]
write.table(values,"./Pearson_values/coding_coding_normalized.txt",sep="\t",quote=F,row.names=F,col.names=F)
# plot densities of the values 
library(ggplot2)
X <- values
data <- data.frame(X)
pdf("./Pearson_graphs/coding_coding_normalized.pdf")
ggplot(data, aes(x = X)) + geom_density() + labs(title = "Neighbors", x="Pearson Correlation", y="Density")
dev.off()
pdf("./Pearson_graphs/coding_coding_normalized_histogram.pdf")
hist(values)
dev.off()

# calc Spearman correlations too
values <- vector(mode="numeric", length=length(coding[,1]))
for (i in 1:length(coding[,1])){
    values[i] <- cor(as.numeric(coding[i,2:100]),as.numeric(coding[i,102:200]),method="spearman")
    if (i%%1000==0) print(i)
}
# remove NAs from values
values <- values[!is.na(values)]
write.table(values,"./Spearman_values/coding_coding_normalized_spearman.txt",sep="\t",quote=F,row.names=F,col.names=F)
pdf("./Spearman_graphs/coding_coding_normalized_spearman.pdf")
ggplot(data, aes(x = X)) + geom_density() + labs(title = "Neighbors", x="Spearman Correlation", y="Density")
dev.off()
pdf("./Spearman_graphs/coding_coding_normalized_histogram_spearman.pdf")
hist(values)
dev.off()

##########################
##### for null model #####
##########################
# null model of random pairs of protein coding genes
neighbors <- c
neighbors <- transform(neighbors[,2], C = sample(neighbors[,4]))
neighbors <- transform(neighbors[,1], C = sample(neighbors[,2]))
neighbors <- neighbors[!duplicated(neighbors[,1]),]
neighbors <- neighbors[!duplicated(neighbors[,2]),]
neighbors <- neighbors[as.character(neighbors[,1])!=as.character(neighbors[,2]),]

# neighbors without duplicate gene_IDs
coding1 <- codingExpr[match(neighbors[,1], codingExpr[,1]),]
# to normalize lncRNA values
pname <- coding1[,1]
pexpr <- coding1[,2:length(coding1)]
plog2 <- log2(pexpr + 1)
psum <- sum(plog2)
pnorm <- plog2 / psum
coding1 <- cbind(pname,pnorm)
# PCG's neighbors
coding2 <- codingExpr[match(neighbors[,2], codingExpr[,1]),]
pname <- coding2[,1]
pexpr <- coding2[,2:length(coding2)]
plog2 <- log2(pexpr + 1)
psum <- sum(plog2)
pnorm <- plog2 / psum
coding2 <- cbind(pname, pnorm)
coding <- cbind(coding1, coding2)
coding <- coding[!is.na(coding[,101]),]
coding <- coding[!is.na(coding[,1]),]
coding <- as.matrix(coding)
test <- cor(as.numeric(coding[,2:100]),as.numeric(coding[,102:200]),method="pearson")
values <- vector(mode="numeric", length=length(coding[,1]))
for (i in 1:length(coding[,1])){
  values[i] <- cor(as.numeric(coding[i,2:100]),as.numeric(coding[i,102:200]),method="pearson")
  if (i%%1000==0) print(i)
}
values <- values[!is.na(values)]
write.table(values,"./Pearson_values/null_model_normalized.txt",sep="\t",quote=F,row.names=F,col.names=F)
library(ggplot2)
X <- values
data <- data.frame(X)
pdf("./Pearson_graphs/null_model_normalized.pdf")
ggplot(data, aes(x = X)) + geom_density() + labs(title = "Neighbors", x="Pearson Correlation", y="Density")
dev.off()
pdf("./Pearson_graphs/null_model_normalized_histogram.pdf")
hist(values)
dev.off()

# calc Spearman correlations too
values <- vector(mode="numeric", length=length(coding[,1]))
for (i in 1:length(coding[,1])){
    values[i] <- cor(as.numeric(coding[i,2:100]),as.numeric(coding[i,102:200]),method="spearman")
    if (i%%1000==0) print(i)
}
values <- values[!is.na(values)]
write.table(values,"./Spearman_values/null_model_normalized_spearman.txt",sep="\t",quote=F,row.names=F,col.names=F)
pdf("./Spearman_graphs/null_model_normalized_spearman.pdf")
ggplot(data, aes(x = X)) + geom_density() + labs(title = "Neighbors", x="Spearman Correlation", y="Density")
dev.off()
pdf("./Spearman_graphs/null_model_normalized_histogram_spearman.pdf")
hist(values)
dev.off()

##################################
## create plot with all 3 lines ##
##################################
# with Pearson values
lncRNA_coding <- read.table("./Pearson_values/lncRNA_coding_normalized.txt",as.is=T,header=F)
coding_coding <- read.table("./Pearson_values/coding_coding_normalized.txt",as.is=T,header=F)
null_model <- read.table("./Pearson_values/null_model_normalized.txt",as.is=T,header=F)

lncRNA_coding$name <- 'lncRNA_coding'
coding_coding$name <- 'coding_coding'
null_model$name <- 'random_coding'

neighbors <- rbind(lncRNA_coding,coding_coding,null_model)
neighbors <- data.frame(neighbors)
pdf("./Pearson_graphs/combined_neighbors_normalized.pdf")
ggplot(neighbors, aes(V1, fill = name)) + geom_density(alpha=0.3) + labs(title = "Neighbors", x="Pearson Correlation", y="Density")
dev.off()

# with Spearman values
lncRNA_coding <- read.table("./Spearman_values/lncRNA_coding_normalized_spearman.txt",as.is=T,header=F)
coding_coding <- read.table("./Spearman_values/coding_coding_normalized_spearman.txt",as.is=T,header=F)
null_model <- read.table("./Spearman_values/null_model_normalized_spearman.txt",as.is=T,header=F)

lncRNA_coding$name <- 'lncRNA_coding'
coding_coding$name <- 'coding_coding'
null_model$name <- 'random_coding'

neighbors <- rbind(lncRNA_coding,coding_coding,null_model)
neighbors <- data.frame(neighbors)
pdf("./Spearman_graphs/combined_neighbors_normalized_spearman.pdf")
ggplot(neighbors, aes(V1, fill = name)) + geom_density(alpha=0.3) + labs(title = "Neighbors", x="Spearman Correlation", y="Density")
dev.off()