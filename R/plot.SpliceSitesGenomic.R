plot.SpliceSitesGenomic <- function(x, col.variant=par("col"), split=FALSE, ...) {
  pos <- x@spsiteIpos
  xlim <- range(0, pos, x@seq.length)
  plot.new()
  plot.window(xlim, c(-1, 1))
  segments(xlim[1], 0, xlim[2], 0, col="grey")
  for (i in seq(1, nrow(pos), length=nrow(pos))) {
    rect(pos[i, 1], -.5, pos[i, 2], +.5)
  }
  col.variant <- rep(col.variant, length=length(x@variants))
  for (i in seq(along=x@variants)) {
    if (split) {
      plot.new()
      plot.window(xlim, c(-1, 1))
      segments(xlim[1], 0, xlim[2], 0, col="grey")
      for (i in seq(1, nrow(pos), length=nrow(pos))) {
        rect(pos[i, 1], -.5, pos[i, 2], +.5)
      }   
    }
    variant <- x@variants[[i]]
    for (j in seq(along=variant)[-1]) {
      x1 <- mean(c(pos[j-1, 1], pos[j-1, 2]))
      x2 <- mean(c(pos[j, 1], pos[j, 2]))
      up <- (-1)^i
      print(c(x1, x2))
                    segments(x1, up*.5, (x1+x2)/2, up*.8, col=col.variant[i])
      segments((x1+x2)/2, up*.8, x2, up*.5, col=col.variant[i])
    }
  }
}
