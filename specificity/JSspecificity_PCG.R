# tissue specificity
file_lncRNA = read.table("anno_lncRNA_all_tissue", as.is=T, sep="\t", header=T)
file <- file_lncRNA
head(file)
dim(file) #3199 rows for now

# normalize V to V' 
V <- function(p) {
  pname <- p[1]
  pexpr <- p[2:length(p)]
  pexpr <- as.numeric(pexpr)
  plog2 <- log2(pexpr + 1)
  psum <- sum(plog2)
  pnorm <- plog2 / psum
  return(pnorm)
}
# H is entropy of discrete prob distribution
H <- function(p) {
  p <- p[p != 0]
  plog <- p*log(p)
  psum <- -sum(plog)
  return(psum)
}
# JS divergence of p1 and p2
JS <- function(p1, p2) {
  pavg <- (p1+p2)/2
  js <- H(pavg) - (H(p1) + H(p2))/2
  return(js)
}
file <- file_lncRNA
# make vector of transcript names
names <- file[,1]
# vector of tissue specificities for each transcript
JS.sp <- vector(mode="numeric", length=length(names))
# calc tissues specificity score for each transcript
for (k in 1:length(names)) {
  p <- file[k,]
  p <- V(p)
  #scores <- vector(mode="numeric", length=length(p))
  o <- matrix(rep(p,each=length(p)),nrow=length(p))
  e <- diag(length(p))
  p_avg <- (o+e)/2
  p_avg.H <- apply(p_avg, 1, H)
  o.H <- apply(o, 1, H)
  e.H <- apply(e, 1, H)
  p_js <- p_avg.H - (o.H + e.H)/2
  p_jsdist <- sqrt(p_js)
  scores <- 1-p_jsdist
  # argmax - maximal tissue specificity score across all n tissues
  if (!all(is.na(scores))) {
    JS.sp[k] = scores[which.max(scores)]
  }
  if (k%%100==0) print(k)
}
js_results <- cbind(names, JS.sp)


