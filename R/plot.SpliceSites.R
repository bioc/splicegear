plot.SpliceSites <- function(x, col.typeI="greenyellow", add=FALSE, ylim=NULL, ...) {

  ## type I splice sites
  if (nrow(x@spsiteIpos) > 0) {
    spliceI.pos <- 1:nrow(x@spsiteIpos)
  } else {
    spliceI.pos <- integer(0)
  }
  
  ## plot the upper part:

  if (! add) {
    xlim <- range(c(1, x@seq.length))
    ylim <- range(0, - nrow(x@spsiteIpos))
    
    plot(x=0, y=0,
         xlab="seq", ylab="splice variants",
         xlim=xlim, ylim=ylim,
         type="n", ...)
  
  } else {
    if (is.null(ylim))
      stop("ylim is missing !")
  }
  
  ## plot type I splice sites
  for (i in spliceI.pos)
    rect(min(x@spsiteIpos[i, ]), ylim[1], max(x@spsiteIpos[i, ]), ylim[2],
         col="yellow", border="transparent" )

  ## plot type II splice sites
  for (i in seq(along=x@spsiteIIpos))
    segments(x@spsiteIIpos[i], ylim[1], x@spsiteIIpos[i], ylim[2], col="red")

  ## plot type III splice sites
  if (nrow(x@spsiteIIIpos) > 0) {
    splice.pos <- 1:nrow(x@spsiteIIIpos)
  } else {
    splice.pos <- integer(0)
  }
  
  for (i in splice.pos) {
    segments(min(x@spsiteIIIpos[i, ]), ylim[1], min(x@spsiteIIIpos[i, ]), ylim[2], col="orange")
    segments(max(x@spsiteIIIpos[i, ]), ylim[1], max(x@spsiteIIIpos[i, ]), ylim[2], col="orange")
  }
  
  ## plot the lower part
  ##FIXME
  ##siteI.ypos <- seq(min(ylim), min(ypos), length=length(spliceI.pos)+1)
  siteI.ypos <- seq(min(ylim), min(c(ylim[2], 0)),
                    length=length(spliceI.pos)+1)
  
  col.typeI <- rep(col.typeI, length=length(spliceI.pos))
  
  for (i in seq(along=spliceI.pos))
    segments(min(x@spsiteIpos[i, ]), siteI.ypos[i],
             max(x@spsiteIpos[i, ]), siteI.ypos[i], col=col.typeI[i])

  ## plot separator
  if (! add)
    abline(h=0, col="grey")
  
  invisible(ylim)
  
}
