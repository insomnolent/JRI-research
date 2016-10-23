setwd("/Users/christinesun/Documents/Research/JRI/gencode/lncRNA_coding")
# read in neighbors
n_lncRNA = read.table("lncRNA_coding_neighbors_revised", as.is=T, header=F)
setwd("/Users/christinesun/Documents/Research/JRI/gencode/coding_coding")
# read in neighbors - has overlap 
n_coding = read.table("coding_coding_neighbors_revised", as.is=T, header=F)

### read in lncRNA:coding neighbors ###
neighbors <- n_lncRNA[n_lncRNA$V4!=n_lncRNA$V11,]
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
# remove rows where transcripts are diff but gene_IDs are the same
neighbors <- dupes
lncRNA_neighbors <- neighbors[as.character(neighbors[,4])!=as.character(neighbors[,12]),]
lncRNA_col3 <- lncRNA_neighbors[,10]-lncRNA_neighbors[,2]
lncRNA_coding <- cbind(lncRNA_neighbors[,4:5],lncRNA_neighbors[12:13],lncRNA_col3)


### read in coding:coding neighbors ###
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
neighbors <- dupes
# get rid of duplicate neighbors (where there's something like gene1:gene2 and also gene2:gene1)
a <- neighbors[neighbors[,1] %in% neighbors[,3] & neighbors[,3] %in% neighbors[,1],]
a <- a[!duplicated(a[,3]),]
b <- neighbors[!(neighbors[,1] %in% neighbors[,3] & neighbors[,3] %in% neighbors[,1]),]
c <- rbind(b,a)
neighbors <- c
coding_neighbors <- neighbors[as.character(neighbors[,4])!=as.character(neighbors[,12]),]
coding_col3 <- coding_neighbors[,10]-coding_neighbors[,2]
coding_coding <- cbind(coding_neighbors[,4:5],coding_neighbors[12:13],coding_col3)


# somehow became factor so need to convert back to character
lncRNA_coding[,5] <- as.numeric(lncRNA_coding[,5])
coding_coding[,5] <- as.numeric(coding_coding[,5])

# take lncRNA:coding and sort into upstream and downstream lncRNAs
# for lncRNA:PCG1 and PCG1:PCG2
lncRNA_PCG1 <- lncRNA_coding[lncRNA_coding[,5]>0,]
lncRNA_left <- lncRNA_PCG1[,1:2]
PCG1_right <- lncRNA_PCG1[,3:4]
PCG2_right <- coding_coding[match(PCG1_right[,2], coding_coding[,2]),]
PCG2_right <- PCG2_right[!is.na(PCG2_right[,1]),]
PCG2_right <- PCG2_right[PCG2_right[,5]>0,]
lncRNA_left <- lncRNA_coding[match(PCG2_right[,2], lncRNA_coding[,4]),]
# for PCG1:lncRNA and PCG2:PCG1
PCG1_lncRNA <- lncRNA_coding[lncRNA_coding[,5]<0,]
lncRNA_right <- PCG1_lncRNA[,1:2]
PCG1_left <- PCG1_lncRNA[,3:4]
PCG2_left <- coding_coding[match(PCG1_left[,2], coding_coding[,2]),]
PCG2_left <- PCG2_left[!is.na(PCG2_left[,1]),]
PCG2_left <- PCG2_left[PCG2_left[,5]<0,]
lncRNA_right <- lncRNA_coding[match(PCG2_left[,2], lncRNA_coding[,4]),]
# combine lncRNA:PCG1 and PCG1:PCG2 columns
lncRNA_PCG <- rbind(lncRNA_left,lncRNA_right)
PCG_PCG <- rbind(PCG2_right, PCG2_left)

# for lncRNA_PCG correlations
coding1 <- noncodingExpr[match(lncRNA_PCG[,2], noncodingExpr[,1]),]
# to normalize lncRNA values
pname <- coding1[,1]
pexpr <- coding1[,2:length(coding1)]
plog2 <- log2(pexpr + 1)
psum <- sum(plog2)
pnorm <- plog2 / psum
coding1 <- cbind(pname,pnorm)
coding2 <- codingExpr[match(lncRNA_PCG[,4], codingExpr[,1]),]
# to normalize PCG values
pnorm <- coding2[,1]
pexpr <- coding2[,2:length(coding2)]
plog2 <- log2(pexpr + 1)
psum <- sum(plog2)
pnorm <- plog2 / psum
coding2 <- cbind(pname, pnorm)

coding <- cbind(coding1, coding2)
coding <- coding[!is.na(coding[,101]),]
coding <- coding[!is.na(coding[,1]),]
coding <- as.matrix(coding)
values <- vector(mode="numeric", length=length(coding[,1]))
for (i in 1:length(coding[,1])){
  values[i] <- cor(as.numeric(coding[i,2:100]),as.numeric(coding[i,102:200]),method="pearson")
  if (i%%1000==0) print(i)
}
lncRNA_PCG_values <- values
length(lncRNA_PCG_values)

# for PCG_PCG correlations
coding1 <- codingExpr[match(PCG_PCG[,2], codingExpr[,1]),]
# to normalize values
pname <- coding1[,1]
pexpr <- coding1[,2:length(coding1)]
plog2 <- log2(pexpr + 1)
psum <- sum(plog2)
pnorm <- plog2 / psum
coding1 <- cbind(pname,pnorm)
coding2 <- codingExpr[match(PCG_PCG[,4], codingExpr[,1]),]
# to normalize values
pnorm <- coding2[,1]
pexpr <- coding2[,2:length(coding2)]
plog2 <- log2(pexpr + 1)
psum <- sum(plog2)
pnorm <- plog2 / psum
coding2 <- cbind(pname, pnorm)
coding <- cbind(coding1, coding2)
coding <- coding[!is.na(coding[,101]),]
coding <- coding[!is.na(coding[,1]),]
coding <- as.matrix(coding)
values <- vector(mode="numeric", length=length(coding[,1]))
for (i in 1:length(coding[,1])){
  values[i] <- cor(as.numeric(coding[i,2:100]),as.numeric(coding[i,102:200]),method="pearson")
  if (i%%1000==0) print(i)
}
PCG_PCG_values <- values
length(PCG_PCG_values)

a <- lncRNA_PCG_values[!is.na(PCG_PCG_values)]
b <- PCG_PCG_values[!is.na(lncRNA_PCG_values)]

lncRNA_PCG_values <- a[!is.na(a)]
PCG_PCG_values <- b[!is.na(b)]

# paired t-test
t.test(lncRNA_PCG_values,PCG_PCG_values,paired=TRUE)

lncRNA_coding <- data.frame(lncRNA_PCG_values)
coding_coding <- data.frame(PCG_PCG_values)
lncRNA_coding$name <- 'lncRNA_coding'
coding_coding$name <- 'coding_coding'

colnames(lncRNA_coding) <- c("value", "name")
colnames(coding_coding) <- c("value", "name")

graph <- rbind(lncRNA_coding,coding_coding)
graph <- data.frame(graph)
ggplot(graph, aes(value, fill = name)) +  scale_x_continuous(limits=c(-1,1)) + theme(panel.background = element_blank())+ geom_density(alpha=0.3) + labs(title = "Neighbors of middle coding transcript", x="Pearson Correlation", y="Density")
