##.initProbes <- function(where=where) {
  
  setClass("Probes",
           representation(pos="matrix", info="data.frame"))
##           where = where)

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
##, where = where)

if( !isGeneric("plot") )
  setGeneric("plot", function(x, y, ...)
             standardGeneric("plot"))
##, where=where)
  
  setMethod("plot",
            signature(x="Probes", y="missing"),
            function(x, ..., probepos.yscale=NULL){ ##, fig.yratio=c(2,1)) {
              plot.Probes(x, ..., probepos.yscale=NULL) ##, fig.yratio=c(2,1))
              })
##            where=where)

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
##, where = where)
  
#   if( !isGeneric("getPos") )
#     setGeneric("getPos", function(object)
#                standardGeneric("getPos"), where=where)
 
#   setMethod("getPos", signature(object="Probes"),
#             function(object, what=c("begin", "end")){
#               what <- match.arg(what)
#               object@info[what][[1]]
#             }, where = where)

##}
