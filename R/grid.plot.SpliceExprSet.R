grid.plot.SpliceExprSet <- function(x, probes.opt = list(),
                               expr.opt = list(col=NA, lty = 1:6),
                               fig.xratio=c(2,1), fig.yratio=c(2,1),
                               probepos.yscale=NULL, ylim = NULL, ...) {
  
  spSites <- x@spliceSites
  probes <- x@probes
  eset <- x@eset

  if (! all(is.list(probes.opt) , is.list(expr.opt)))
    stop("probes.opt and expr.opt should be lists !")
    
  if (all(is.na(expr.opt$col))) {
    expr.opt$col <- rainbow(ncol(exprs(eset)))
  }
  
  grid.newpage()
  ##top.vp <- top.vp <- viewport(y = 0, height = unit(1, "npc") - unit(1.5, 
  ##      "lines"), just = c("centre", "bottom"))
  top.lt <- grid.layout(2, 2, widths = 0.9 * c(0.5),
                        heights = 0.9 * 1, ##c(0.5),
                        default.units = "npc",
                        respect = matrix(c(1, 1), 2, 2))
  temp.vp <- viewport(layout = top.lt)
  push.viewport(temp.vp)

  panel.vp <- viewport(layout.pos.row = 2, layout.pos.col = 1)  
  grid.plot(spSites, vp=panel.vp, add=TRUE)
  
  panel.vp <- viewport(layout.pos.row = 1, layout.pos.col = 1)
  ylim <- do.call("grid.plot", c(list(probes, xlim=c(0, spSites@seq.length), vp=panel.vp, add=TRUE),
                                 probes.opt))
  
  
  if (is.null(probepos.yscale)) {
    ypos <- 1:nrow(probes@pos)
  } else {
    ypos <- probepos.yscale
  }

  ##ylim <- range(ypos)
    
  if (is.null(ylim))
    ylim <- m.ylim
  
#     if (nrow(probes@pos) <= 1)
#     ylim <- c(0,1)
#   else
#     ylim <- c(0,nrow(probes@pos))

  
  if (nrow(exprs(eset)) != nrow(probes@pos))
    stop("length mismatch between number of probes and number of expression values")
  
  opar.mar <- par()$mar
  on.exit(par(mar=opar.mar))
  npar.mar <- opar.mar
  npar.mar[2] <- 1
  par(mar=opar.mar)

  ylim <- range(ylim, ylim[1] - 1/4 * (ylim[2] - ylim[1]))
  xlim <- c(0, spSites@seq.length)

  ##scale.x <- grid.make.numeric2npc(xlim)
  ##scale.y <- grid.make.numeric2npc(ylim)

  gp <- grid.expand.gp(nrow(exprs(eset)), parlist=expr.opt)

  panel.vp <- viewport(layout.pos.row = 1, layout.pos.col = 2)
  push.viewport(panel.vp)
  vp <- viewport(xscale = xlim, yscale = ylim)
  for (i in seq(1, nrow(exprs(eset)), length=nrow(exprs(eset))))
    grid.lines(exprs(eset)[i, ], ypos,
               gp = gp[[i]], vp = panel.vp)
  
  return()
  
  do.call("matplot", c(list(exprs(eset), matrix(ypos, ncol=1),
                            ylim=ylim,
                            xlab="expression", ylab="probes",
                            type="l"), expr.opt))
  ##overlay typeI
                                        #for (i in 1:nrow(spSites@spsiteIpos))
                                        #  rect(ylim[1], min(spSites@spsiteIpos[i, ]), ylim[2], max(spSites@spsiteIpos[i, ]),
                                        #       col="yellow", border="transparent" )
  ##overlay typeII
                                        #for (i in seq(along=spSites@spsiteIIpos))
                                        #  ##segments(x@spsiteIIpos[i], ylim[1], x@spsiteIIpos[i], ylim[2], col="red")
                                        #  segments(ylim[1], spSites@spsiteIIpos[i], ylim[2], spSites@spsiteIIpos[i], col="red")
  
  ##boxplot(exprs(eset), horizontal=TRUE, add=TRUE)
  
}
