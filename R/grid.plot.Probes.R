
grid.plot.Probes <- function(x,
                             col="black",
                             add=FALSE,
                             probepos.yscale=NULL, xlim=NULL, vp = NULL, ...) {

  if (is.null(probepos.yscale)) {
    if (nrow(x@pos) <= 1)
      ypos <- c(0,1)
    else
      ypos <- 1:nrow(x@pos)
  } else {
    ypos <- probepos.yscale
  }

  ylim <- range(ypos)
  
  if (! add) {
    grid.newpage()

    ##grid.rect(gp=gpar(fill="grey"))
    if (is.null(xlim))
      xlim <- range(x@pos)
    
    ## plot separator
    ##abline(h=0, col="grey")
    vp <- viewport(xscale = xlim, yscale = ylim, w=0.9, h=0.9)
    
  } else {
    pushViewport(vp)
    on.exit(popViewport())
    vp <- viewport(xscale = xlim, yscale = ylim)
  }
  
  
  grid.xaxis(vp = vp)

  if (nrow(x@pos) > 0) {
    col <- rep(col, length=length(ypos))
    for (i in seq(along=ypos)) {
      grid.segments(x@pos[i, 1], ypos[i],
                    x@pos[i, 2], ypos[i],
                    default.units = "native",
                    gp = gpar(col=col[i]),
                    vp = vp)
##      segments(x@pos[i, 1], ypos[i], x@pos[i, 2], ypos[i], col=col[i])
    }
  }

  invisible(ylim)
  
}
