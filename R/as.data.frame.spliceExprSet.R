as.data.frame.SpliceExprSet <- function(x, row.names=NA, optional=NA) {
  
  if (! inherits(x, "SpliceExprSet"))
    stop("the argument should inherit class 'SpliceExprSet'")

  spSites <- x@spliceSites
  eset <- x@eset
  probes <- x@probes
  
  nc.eset <- ncol(exprs(eset))
  nr.eset <- nrow(exprs(eset))
  
  ## build the probe position column
  i.probes <- seq(1, nr.eset)
  i.probes.expand <- rep(i.probes, nc.eset)
  
  ##ppos1 <- rep(probes@pos[, 1], rep(nc.eset, length(probes@pos[, 1])))
  ##ppos2 <- rep(probes@pos[, 2], rep(nc.eset, length(probes@pos[, 2])))
  
  ## build the 'is in type I' column
  probeOnSpSite <- isProbeOnSpliceSite(probes, spSites)

  ##probeOnSpSite <- lapply(probeOnSpSite,
  ##                        function(x) rep(x, rep(nc.eset, length(x)))
  ##                        )
  
  ## build the probe level intensities columns
  rv.eset <- as.data.frame.exprSet(eset)
  ## expand the covariate in phenoData
  
  ## build the covariate information
  ##FIXME: to be implemented...
  
  # rv <- do.call("data.frame", c(list(ppos1=ppos1, ppos2=ppos2,
#                                      isintypeI=probeOnSpSite$isintypeI,
#                                      isintypeII=probeOnSpSite$isintypeII,
#                                      exprs=pexp),
#                                      )

  rv <- do.call("data.frame", c(list(begin = probes@pos[, 1][i.probes.expand],
                                     end = probes@pos[, 2][i.probes.expand],
                                     isintypeI = probeOnSpSite$isintypeI[i.probes.expand],
                                     isintypeII = probeOnSpSite$isintypeII[i.probes.expand]),
                                lapply(rv.eset, function(x) {return(x)})
                                )
                )
                
  return(rv)
}
