\name{plot.SpliceExprSet}
\alias{plot.SpliceExprSet}
\title{ plot a SpliceExprSet }
\description{
  Plot a object of class \code{SpliceExprSet}
}
\usage{
plot.SpliceExprSet(x, ..., probes.opt = list(), expr.opt = list(col = NA, lty = 1:6), fig.xratio = c(2, 1), fig.yratio = c(2, 1), probepos.yscale = NULL)
}x
\arguments{
  x\item{x}{ a \code{\link{SpliceExprSet-class}}}
  \item{\dots}{ optional parameters to be passed to the function \code{plot} }
  \item{probes.opt}{ optional parameters to be passed for the plotting of
  the \code{\link{Probes-class}}}
  \item{expr.opt}{ optional parameters to be passed for the plotting of
  the \code{\link{exprSet-class}}}
  \item{fig.xratio}{ ratio between the left and right parts of the plot }
  \item{fig.yratio}{ ratio between the upper and lower parts of the plot }
  \item{probepos.yscale}{ enforce `y' positions for the probes. }
}
\details{

}
\value{
  function used for its side-effect(s).
}

\author{ laurent }

\seealso{\code{\link{SpliceExprSet-class}} }
\examples{
data(spliceset.Hs.44)

speset <- spliceset.Hs.44

levels(pData(speset@eset)$Material)
## Liver, Mix and SNB19
cl.mat <- c("red", "yellow","blue")[as.integer(pData(speset@eset)$Material)]
## colored in red, yellow and blue respectively
plot(speset, expr.opt=list(col=cl.mat, log="x"))

begin.pos <- speset@probes@pos[,1]
plot(speset, expr.opt=list(col=cl.mat), probespos.yscale=begin.pos)

}
\keyword{ hplot }
