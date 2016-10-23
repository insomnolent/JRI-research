# tissue specificity
file = read.table("/Users/christinesun/Documents/Research/JRI/JSdata.txt", as.is=T, sep="\t", header=T)
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

# make vector of transcript names
names <- file[,1]
# vector of tissue specificities for each transcript
JS.sp <- vector(mode="numeric", length=length(names))
# calc tissues specificity score for each transcript
for (k in 1:length(names)) {
  p <- file[k,]
  p <- V(p)
  scores <- vector(mode="numeric", length=length(p))
  for(i in 1:length(p)) {
    # e.t is a predefined expr pattern where a transcript is expressed in only one tissue
    e <- vector(mode="numeric", length=length(p))
    e[i] = 1
    p.js <- JS(p,e)
    p.jsdist <- sqrt(p.js)
    scores[i] = 1 - p.jsdist
  }
  # argmax - maximal tissue specificity score across all n tissues
  if (!all(is.na(scores))) {
    JS.sp[k] = scores[which.max(scores)]
  }
  # if transcript has no expression in any tissue
  # JS.sp[k] = NA
}
js_results <- cbind(names, JS.sp)
lincRNA <- 