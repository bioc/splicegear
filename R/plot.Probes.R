
plot.Probes <- function(x,
                        col="black",
                        add=FALSE,
                        probepos.yscale=NULL, xlim=NULL, ...) {

  if (is.null(probepos.yscale)) {
    if (nrow(x@pos) <= 1)
      ypos <- c(0,1)
    else
      ypos <- 1:nrow(x@pos)
  } else {
    ypos <- probepos.yscale
  }

  ylim <- numeric(0)
  
  if (! add) {
    if (is.null(xlim)) {
      if (nrow(x@pos) == 0)
        xlim <- c(0,1)
      else
        xlim <- range(x@pos)
    }
    ## plot the upper part:
    
    ylim <- range(ypos, 0)
    
    plot(x=0, y=0,
         xlab="seq", ylab="probes",
         xlim=xlim, ylim=ylim,
         type="n", ...)
    
    ## plot separator
    abline(h=0, col="grey")

  }
  
  if (nrow(x@pos) > 0) {
    col <- rep(col, length=length(ypos))
    for (i in seq(along=ypos)) {
      segments(x@pos[i, 1], ypos[i], x@pos[i, 2], ypos[i], col=col[i])
    }
  }

  if (identical(ylim, numeric(0)))
    ylim <- ypos
  else
    ylim <- range(ypos)
  
  invisible(ylim)
  
}
