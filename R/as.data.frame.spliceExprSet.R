as.data.frame.SpliceExprSet <- function(x, row.names=NA, optional=NA, ...) {

    if (! is(x, "SpliceExprSet"))
        stop("the argument should inherit class 'SpliceExprSet'")

  spSites <- x@spliceSites
  eset <- x@eset
  probes <- x@probes

  nc.eset <- ncol(exprs(eset))
  nr.eset <- nrow(exprs(eset))

  ## sanity check
  ## (because slots of objects can be tweaked by the user)
  if (nr.eset != nrow(probes@pos)) {
    stop("inconsistency beetween slots 'eset' and 'probes'")
  }

  ## build the probe position column
  i.probes <- seq(1, nr.eset)
  i.probes.expand <- rep(i.probes, nc.eset)

  ## build the 'is in type x' columns
  probeOnSpSite <- isProbeOnSpliceSite(probes, spSites)

  ## build the probe level intensities columns
  rv.eset <- local({
    X <- exprs(eset)
    cbind(as.data.frame(lapply(pData(eset), rep, each=nrow(X))),
                        exprs=c(X))
        })

  rv <- do.call(data.frame, c(list(begin = probes@pos[, 1][i.probes.expand],
                                     end = probes@pos[, 2][i.probes.expand],
                                     isintypeI = probeOnSpSite$isintypeI[i.probes.expand],
                                     isintypeII = probeOnSpSite$isintypeII[i.probes.expand]),
                                lapply(rv.eset, function(x) return(x)),
                                lapply(probes@info, function(x, i) return(x[i]), i.probes.expand)
                                )
                )

  return(rv)
}
