plot.SpliceExprSet <- function(x, ..., probes.opt = list(), expr.opt = list(col=NA, lty = 1:6),
                               fig.xratio=c(2,1), fig.yratio=c(2,1),
                               probepos.yscale=NULL) {
  
  ##               plotOneProbeWithColor <- function(posx, y, exprs) {
  ##                 segments(posx[1], y[1], posx[2], y[1], ...)                
  ##               }
  spSites <- x@spliceSites
  probes <- x@probes
  eset <- x@eset

  if (! all(is.list(probes.opt) , is.list(expr.opt)))
    stop("probes.opt and expr.opt should be lists !")
    
  if (all(is.na(expr.opt$col))) {
    expr.opt$col <- rainbow(ncol(exprs(eset)))
  }
  
  layout(matrix(c(1,2), 1, 2), width=fig.xratio, height=fig.yratio) #fig.yratio useless for now
  
  if (is.null(probepos.yscale)) {
    ypos <- 1:nrow(probes@pos)
  } else {
    ypos <- probepos.yscale
  }

  ##ylim <- range(ypos)

  plotSplSites <- getMethod("plot", c("Probes", "SpliceSites"))
  plotSplSites(probes, spSites, ..., fig.yratio=fig.yratio, probepos.yscale=probepos.yscale)
  ##ylim <- plot(spSites, ..., fig.yratio=fig.yratio, probepos.yscale=probepos.yscale)
  if (nrow(exprs(eset)) != nrow(probes@pos))
    stop("length mismatch between number of probes and number of expression values")
  
  opar.mar <- par()$mar
  on.exit(par(mar=opar.mar))
  npar.mar <- opar.mar
  npar.mar[2] <- 1
  par(mar=opar.mar)
  if (nrow(probes@pos) <= 1)
    ylim <- c(0,1)
  else
    ylim <- c(0,nrow(probes@pos))
  ylim <- range(ylim, ylim[1] - 1/4 * (ylim[2] - ylim[1]))
  xlim <- c(0, spSites@seq.length)

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
