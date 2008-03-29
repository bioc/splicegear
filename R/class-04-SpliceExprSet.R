## copyright. Laurent 2003
## under the L-GPL license

## classe to describe splice variants w/ expression values
##
require(Biobase, quietly=TRUE) || stop("Could not load the package 'Biobase'")

setClass("SpliceExprSet",
         representation(spliceSites="SpliceSites",
                        probes = "Probes",
                        eset="ExpressionSet")) ## exprs: one row per probe,
## one column per experiment
## AnnotatedDataFrame: covariate info.

## -- accessors --
if (is.null(getGeneric("spliceSites")))
  setGeneric("spliceSites", function(object)
               standardGeneric("spliceSites"))

setMethod("spliceSites", signature(object = "SpliceExprSet"),
          function(object) {
            object@spliceSites
          })


setMethod("exprs", signature(object="SpliceExprSet"),
          function(object) exprs(object@eset))

setReplaceMethod("exprs", "SpliceExprSet",
                 function(object, value) {
                   exprs(object@eset) <- value
                   if (nrow(exprs(object)) != nrow(object@probes@pos))
                     stop("mismatch between the number of probes and the size of the expression matrix.")
                   return(object)
                 })

##setMethod("probepos", signature(object="SpliceExprSet"),
##          function(object) { return(object@probepos) }, where=where)

##setMethod("variant", signature(object="SpliceExprSet"),
##          function(object) { return(object@variant) }, where=where)

## -- other methods --
setMethod("show", signature(object = "SpliceExprSet"),
          function(object) {
            cat("Alternative splicing expression set (SpliceExprSet):\n")
            cat("\t", ncol(exprs(object@eset)), " chip(s)\n", sep="")
            cat("\t", nrow(exprs(object)), " probe(s) on the sequence\n", sep="")
            show(spliceSites(object))
          })


setMethod("plot", signature(x = "SpliceExprSet", y = "missing"),
          function(x, probes.opt = list(), expr.opt = list(col=NA, lty = 1:6),
                   fig.xratio=c(2,1), fig.yratio=c(2,1),
                   probepos.yscale=NULL, ...) {
            plot.SpliceExprSet(x, probes.opt = probes.opt,
                               expr.opt = expr.opt,
                               fig.xratio=fig.xratio, fig.yratio=fig.yratio,
                               probepos.yscale=probepos.yscale, ...)
          })

if( !isGeneric("grid.plot") )
  setGeneric("grid.plot", function(x, y, ...)
             standardGeneric("grid.plot"))

setMethod("grid.plot",
          signature(x="SpliceExprSet", y="missing"),
          function(x, ...) {
            grid.plot.SpliceExprSet(x, ...)
          })


##setMethod("sort", signature(x = "SpliceExprSet"),
sort.SpliceExprSet <- function(x, decreasing, fun=function(x) order(x@probes@pos[, 1]), reverse=FALSE, ...)
{
  
  o <- fun(x)
  
  if (reverse) {
    o <- rev(o)
  }
  
  spSites <- x@spliceSites
  eset <- x@eset
  probes <- x@probes
  probes@pos <- probes@pos[o, , drop=FALSE]
  exprs(eset) <- exprs(eset)[o, , drop=FALSE]
  ##if (! is.null(spSites@variant))
  ##  spSites@variant <- spSites@variant[o]

  x@spliceSites <- spSites
  x@eset <- eset
  x@probes <- probes
  
  return(x)
}

##  setMethod("isProbeOnSpliceSite", "SpliceExprSet",
isSpliceSiteOnProbe <- function(spSites, probes) {
  stop("Not implemented (yet)")
}
##, where=where)

isProbeOnSpliceSite <- function(probes, spSites) {
  
  ## ensure that pos1 & pos2 are correctly ordered
  r.ppos <- apply(probes@pos, 1, range)
  
  isWithin <- function(x,y) {
    ## x: pos1 and pos2 for a probe
    ## y: matrix of two rows and n columns (pos1 and pos2 for n type-II splice sites)
    ## check if (pos1 of x is greater the pos1 of y and lower than the value pos2 of y
    ## for any column (i.e. for any type-II splice site) of y.
    ## In other words, it returns TRUE if the probe x has at least a partial overlap with
    ## one of the type-II site, FALSE otherwise
    any(x[1] >= y[1, ] & x[1] <= y[2, ]) | any(x[2] <= y[2, ] & x[2] >= y[1, ])
  }

  if(length(spSites@spsiteIpos) > 0)
    isintypeI <- apply(r.ppos, 2, isWithin,
                       apply(spSites@spsiteIpos, 1, range))
  else
    isintypeI <- rep(FALSE, nrow(probes@pos))
  
  hasSite <- function(x, y) {
    any(x[1] >= y & x[2] <= y)
  }

  if(length(spSites@spsiteIIpos) > 0)
    isintypeII <- apply(r.ppos, 2, hasSite,
                        spSites@spsiteIIpos)
  else
    isintypeII <- rep(FALSE, nrow(probes@pos))
  
  r <- list(isintypeI=isintypeI, isintypeII=isintypeII)
  
  return(r)
}

## HOWTO:
## - ProbeSet to SpliceExprSet:
##            get the pm from the probe set and put in slot exprs
##
