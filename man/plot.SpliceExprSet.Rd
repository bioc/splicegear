\name{plot.SpliceExprSet}
\alias{plot.SpliceExprSet}
\title{ plot a SpliceExprSet }
\description{
  Plot a object of class \code{SpliceExprSet}
}
\usage{
\method{plot}{SpliceExprSet}(x,
                   probes.opt = list(), expr.opt = list(col = NA, lty = 1:6),
                   fig.xratio = c(2, 1), fig.yratio = c(2, 1),
                   probepos.yscale = NULL, ylim,
                   ...)
}
\arguments{
  \item{x}{ a \code{\link{SpliceExprSet-class}}}
  \item{probes.opt}{ optional parameters to be passed for the plotting of
    the \code{\link{Probes-class}} }
  \item{expr.opt}{ optional parameters to be passed for the plotting of
    the \code{\link[Biobase]{ExpressionSet-class}} }
  \item{fig.xratio}{ ratio between the left and right parts of the plot }
  \item{fig.yratio}{ ratio between the upper and lower parts of the plot }
  \item{probepos.yscale}{ enforce `y' positions for the probes. }
  \item{ylim}{range for the y-axis}
  \item{\dots}{ optional parameters to be passed to the function \code{plot} }
}
\details{
  The argument \code{probepos.yscale} can be used to scale probes
  according to their position on the reference sequence, as shown
  in the last example below.
}
\value{
  function used for its side-effect(s).
}

\author{ laurent }

\seealso{\code{\link{SpliceExprSet-class}} }
\examples{
data(spliceset)


levels(pData(spliceset@eset)$Material)
## Liver, Mix and SNB19
cl.mat <- c("red", "yellow","blue")[as.integer(pData(spliceset@eset)$Material)]
## colored in red, yellow and blue respectively
plot(spliceset, expr.opt = list(col = cl.mat, log = "x"))

## sort
spliceset <- sort.SpliceExprSet(spliceset)
begin.pos <- spliceset@probes@pos[, 1]
plot(spliceset, expr.opt = list(col=cl.mat), probepos.yscale = begin.pos)

}
\keyword{ hplot }
