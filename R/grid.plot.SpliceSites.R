grid.plot.SpliceSites <- function(x, col.typeI="orange", col.typeI.window="yellow", col.typeII="red", add=FALSE, ylim=NULL, vp = NULL, ...) {

  ## type I splice sites
  if (nrow(x@spsiteIpos) > 0) {
    spliceI.pos <- 1:nrow(x@spsiteIpos)
  } else {
    spliceI.pos <- integer(0)
  }
  
  ## plot the upper part:

  xlim <- range(c(1, x@seq.length))
  ylim <- range(0, - nrow(x@spsiteIpos))

  if (! add) {
    grid.newpage()
    vp <- viewport(xscale = xlim, yscale = ylim, w=0.9, h=0.9)
    grid.xaxis(vp = vp, main=FALSE)
  } else {
    pushViewport(vp)
    on.exit(popViewport())
    vp <- viewport(xscale = xlim, yscale = ylim)
  }
  
  
  scale.x <- grid.make.numeric2npc(xlim=xlim)
  scale.y <- grid.make.numeric2npc(xlim=ylim)
  ## plot type I splice sites
  for (i in spliceI.pos)
    grid.rect(scale.x(min(x@spsiteIpos[i, ])), ##scale.y(ylim[1]),
              width = scale.x(diff(range(x@spsiteIpos[i, ]))),
              height = 0.99,
              just ="left",
              gp = gpar(fill=col.typeI.window, col="transparent"),
              vp = vp)
  
  ## plot type II splice sites
  for (i in seq(along=x@spsiteIIpos))
    grid.segments(scale.x(x@spsiteIIpos[i]), scale.y(ylim[1]),
                  scale.x(x@spsiteIIpos[i]), scale.y(ylim[2]),
                  gp = gpar(col=col.typeII),
                  vp = vp)
  
#   ## plot type III splice sites
#   if (nrow(x@spsiteIIIpos) > 0) {
#     splice.pos <- 1:nrow(x@spsiteIIIpos)
#   } else {
#     splice.pos <- integer(0)
#   }
  
#   for (i in splice.pos) {
#     grid.segments(scale.x(min(x@spsiteIIIpos[i, ])), scale.y(ylim[1]),
#                   scale.x(min(x@spsiteIIIpos[i, ])), scale.y(ylim[2]),
#                   gp = gpar(col="orange"),
#                   vp = vp)
#     grid.segments(scale.x(max(x@spsiteIIIpos[i, ])), scale.y(ylim[1]),
#                   scale.x(max(x@spsiteIIIpos[i, ])), scale.y(ylim[2]),
#                   gp = gpar(col="orange"),
#                   vp = vp)
#   }
  
  ## plot the lower part
  ##FIXME
  ##siteI.ypos <- seq(min(ylim), min(ypos), length=length(spliceI.pos)+1)
  siteI.ypos <- seq(min(ylim), min(c(ylim[2], 0)),
                    length=length(spliceI.pos)+1)
  
  col.typeI <- rep(col.typeI, length=length(spliceI.pos))
  
  for (i in seq(along=spliceI.pos))
    grid.segments(scale.x(min(x@spsiteIpos[i, ])), scale.y(siteI.ypos[i]),
                  scale.x(max(x@spsiteIpos[i, ])), scale.y(siteI.ypos[i]),
                  gp = gpar(col=col.typeI[i]),
                  vp = vp)

  ## plot separator
  ##if (! add)
  ##  abline(h=0, col="grey")
  
  invisible(ylim)
  
}
