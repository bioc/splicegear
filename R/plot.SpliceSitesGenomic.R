plot.SpliceSitesGenomic <- function(x, col.variant=par("col"), col.exon="white",
                                    split=FALSE, main=names(x@variants), ...) {
  makeBackground <- function(xlim, main) {
    plot.new()
    plot.window(xlim, c(-1, 1))
    title(main=main)
    segments(xlim[1], 0, xlim[2], 0, col="grey")
    for (i in seq(1, nrow(pos), length=nrow(pos))) {
      rect(pos[i, 1], -.5, pos[i, 2], +.5, col=col.exon[i])
    }
  }
  pos <- x@spsiteIpos
  col.variant <- rep(col.variant, length=length(x@variants))
  col.exon <- rep(col.exon, length=nrow(pos))
  xlim <- range(0, pos, x@seq.length)
  if (split)
    main <- paste("variant ", rep(main, length=length(x@variants)))
  else
    makeBackground(xlim, main)
  for (i in seq(along=x@variants)) {
    if (split) {
      makeBackground(xlim, main=main[i])
    }
    variant <- x@variants[[i]]
    for (j in seq(along=variant)[-1]) {
      x1 <- mean(c(pos[variant[j-1], 1], pos[variant[j-1], 2]))
      x2 <- mean(c(pos[variant[j], 1], pos[variant[j], 2]))
      up <- (-1)^i
      segments(x1, up*.5, (x1+x2)/2, up*.8, col=col.variant[i])
      segments((x1+x2)/2, up*.8, x2, up*.5, col=col.variant[i])
    }
  }
}
