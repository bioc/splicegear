## Laurent 2003
## Implemented from Genome Research, Hu et. al. , p.1244


## x <- matrix(1:20, 5, 4); tissue <- factor(c("a", "a", "b", "b"))

## SPLICE algorithm:
getRelSignStrength <- function(x, tissue=as.factor(1:ncol(x)), fun=mean, nipt=30, nitt=30, ...) {
  ## x: matrix. One probe per line, one column per chip
  ## tissue: a factor for the 'tissue' covariate
  ## returns the 'RSS'
  
  ## cutoff values
  x[ x < 20 ] <- 20
  x[ x > 5000 ] <- 5000

  ## mean (or 'fun') of all probe pairs in tissue X 
  avgDix.perx <- unlist(tapply(seq(along=tissue), tissue, function(y) fun(x[ , y])))
  avgDix.perx.indices <- tapply(seq(along=tissue), tissue, function(y) y)
  avgDix <- matrix(NaN, nc=ncol(x), nr=nrow(x))
  for (i in seq(along=avgDix.perx.indices))
    avgDix[, avgDix.perx.indices[[i]]] <- avgDix.perx[i]

  ## mean or 'fun' of a particular probe pair across different tissues
  ##avgDi <- t(apply(x, 1, function(y) tapply(y, tissue, fun)), ...)
  avgDi <- apply(x, 1, fun, ...)
  
  relsignstr <- x / avgDix[, as.integer(tissue)]

  ## non-informative probe threshold
  nip <- avgDi < nipt
  attr(relsignstr, "nip") <- nip

  ## non-informative tissue threshold
  nit <- avgDix < nitt
  attr(relsignstr, "nit") <- nit

  return(relsignstr)
}

getFinalRatio <- function(x, tissue=as.factor(1:ncol(x)), fun=mean, ...) {
  x <- getRelSignStrength(x, tissue=tissue, fun=fun, ...)
  ex <- sapply(tissue, function(x) tissue != x, simplify=FALSE)
  avgrss <- t(apply(x, 1, function(y) unlist(lapply(ex, function(z) fun(y[z])))))
  fr <- log(x / avgrss)

  return(fr)
}

