## This class could be useful in the pack matchprobes too...

setClass("Probes",
         representation(pos="matrix", info="data.frame"))


setMethod("initialize", "Probes",
          function(.Object, pos, info=NULL) {
            .Object@pos <- pos
            if (! is.null(info)) {
              if (nrow(pos) != nrow(info))
                stop("length mismatch between 'pos' and 'info'.")
              .Object@info <- info
            }
            return(.Object)
          })

if( !isGeneric("plot") )
  setGeneric("plot", function(x, y, ...)
             standardGeneric("plot"))

setMethod("plot",
          signature(x="Probes", y="missing"),
          function(x, ..., probepos.yscale=NULL){ ##, fig.yratio=c(2,1)) {
            plot.Probes(x, ..., probepos.yscale=NULL) ##, fig.yratio=c(2,1))
          })

setMethod("plot",
          signature(x="Probes", y="SpliceSites"),
          function(x, y, probes.opt=list(), spsites.opt=list(), ...) {
            if (nrow(x@pos) <= 1)
              ylim <- c(0,1)
            else
              ylim <- c(0,nrow(x@pos))
            ylim <- range(ylim, ylim[1] - 1/4 * (ylim[2] - ylim[1]))
            xlim <- c(0, y@seq.length)
            
            plot(x=0, y=0,
                 xlab="seq", ylab="splice variants / probes",
                 xlim=xlim, ylim=ylim,
                 type="n", ...)
            
            do.call("plot.SpliceSites", c(list(y, add=TRUE, ylim=ylim), spsites.opt))
            do.call("plot.Probes", c(list(x, xlim=c(0, y@seq.length), add=TRUE), probes.opt))
            
            abline(h=0, col="grey")
              
          })


matchprobes2Probes <- function(mpo, probes.length, ref.seq=NULL, ref.seq.length=NULL, names=NULL) {
  if (! identical(names(mpo), c("match", "pos")))
    stop("Expected a list with names 'match' and 'probes'\n(as returned by the package 'matchprobes').")

  n.seq <- length(mpo$pos)
  
  p.list <- vector("list", length=n.seq)

  for (i in seq(along=mpo$pos)) {
    info <- data.frame(probe.index = mpo$match[[i]])
    
    p.list[[i]] <- new("Probes", pos=cbind(mpo$pos[[i]],
                                   mpo$pos[[i]] + probes.length),
                       info=info)
    
  }
  
  return(p.list)
}


#   if( !isGeneric("getPos") )
#     setGeneric("getPos", function(object)
#                standardGeneric("getPos"), where=where)
 
#   setMethod("getPos", signature(object="Probes"),
#             function(object, what=c("begin", "end")){
#               what <- match.arg(what)
#               object@info[what][[1]]
#             }, where = where)

##}
