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

setMethod("show", "Probes",
          function(object) {
            cat("Probes object:\n")
            cat(" ", nrow(object@pos), "probe(s)\n")
          })

if( !isGeneric("plot") )
  setGeneric("plot", function(x, y, ...)
             standardGeneric("plot"))

setMethod("plot",
          signature(x="Probes", y="missing"),
          function(x, probepos.yscale=NULL, ...){ ##, fig.yratio=c(2,1)) {
            plot.Probes(x, probepos.yscale=NULL, ...) ##, fig.yratio=c(2,1))
          })

setMethod("plot",
          signature(x="Probes", y="SpliceSites"),
          function(x, y, probes.opt=list(), spsites.opt=list(), fig.yratio = c(2,1), probepos.yscale = NULL, ...) {
            
            if (is.null(probepos.yscale)) {
              if (nrow(x@pos) <= 1)
                ypos <- c(0,1)
              else
                ypos <- 1:nrow(x@pos)
            } else {
              ypos <- probepos.yscale
            }
            
            ylim <- range(ypos, 0)
            
            ylim <- range(ylim, ylim[1] - 1/4 * (ylim[2] - ylim[1]))
            xlim <- c(0, y@seq.length)
            
            plot(x=0, y=0,
                 xlab="seq", ylab="splice variants / probes",
                 xlim=xlim, ylim=ylim,
                 type="n", ...)
            
            do.call("plot.SpliceSites", c(list(y, add=TRUE, ylim=ylim), spsites.opt))
            p.ylim <- do.call("plot.Probes",
                              c(list(x, xlim=c(0, y@seq.length), add=TRUE, probepos.yscale = probepos.yscale),
                                probes.opt))
            
            abline(h=0, col="grey")


            invisible(ylim)
            
          }
          )


if( !isGeneric("grid.plot") )
  setGeneric("grid.plot", function(x, y, ...)
             standardGeneric("grid.plot"))

setMethod("grid.plot",
          signature(x="Probes", y="missing"),
          function(x, probepos.yscale=NULL, vp = NULL, ...){ ##, fig.yratio=c(2,1)) {
            grid.plot.Probes(x, probepos.yscale=NULL, vp = vp, ...) ##, fig.yratio=c(2,1))
          })


setMethod("grid.plot",
          signature(x="Probes", y="SpliceSites"),
          function(x, y, probes.opt=list(), spsites.opt=list(), fig.yratio = c(2/3, 1/3),
                   probepos.yscale = NULL, add=FALSE, vp = NULL, ...) {
            
            if (is.null(probepos.yscale)) {
              if (nrow(x@pos) <= 1)
                ypos <- c(0,1)
              else
                ypos <- 1:nrow(x@pos)
            } else {
              ypos <- probepos.yscale
            }
            
            ylim <- range(ypos, 0)
            
            ylim <- range(ylim, ylim[1] - 1/4 * (ylim[2] - ylim[1]))
            xlim <- c(0, y@seq.length)



            if (! add) {
              grid.newpage()
              figscale <- 0.9 
              ##vp <- viewport(xscale = xlim, yscale = ylim, w=0.9, h=0.9)
            } else {
              push.viewport(vp)
              on.exit(pop.viewport())
              figscale <- 1
            }
            
            top.lt <- grid.layout(2, 1, widths = figscale * 1,
                                    heights = figscale * fig.yratio,
                                  default.units = "npc",
                                  respect = matrix(c(1, 1), 2, 1))
            
            temp.vp <- viewport(layout = top.lt)
            push.viewport(temp.vp)

            ##spliceSites
            panel.vp <- viewport(layout.pos.row = 2, layout.pos.col = 1)
            do.call("grid.plot", c(list(y, add=TRUE, ylim=ylim, vp=panel.vp), spsites.opt))
            ##probes
            panel.vp <- viewport(layout.pos.row = 1, layout.pos.col = 1)
            ## trick to have the background:
            spsites.opt.hack <- spsites.opt
            spsites.opt.hack$col.typeI <- 0
            do.call("grid.plot", c(list(y, add=TRUE, ylim=ylim, vp=panel.vp),
                                   spsites.opt.hack))
            p.ylim <- do.call("grid.plot",
                              c(list(x, xlim=c(0, y@seq.length), add=TRUE, vp=panel.vp,
                                     probepos.yscale = probepos.yscale),
                                probes.opt))
            
            ##abline(h=0, col="grey")
            
            invisible(ylim)
            
          }
          )



matchprobes2Probes <- function(mpo, probes.length, names=NULL) {
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
