\name{grid.plot.Probes}
\alias{grid.plot.Probes}
\alias{grid.plot.SpliceSites}
\alias{grid.plot.SpliceExprSet}

\title{ Plot splicegear objects }
\description{
  Plot objects defined in the package splicegear
}
\usage{
grid.plot.Probes(x, col = "black", add = FALSE, probepos.yscale = NULL,
              xlim = NULL, vp = NULL, ...)

grid.plot.SpliceSites(x, col.typeI = "orange", col.typeI.window = "yellow", 
    col.typeII = "red", add = FALSE, ylim = NULL, vp = NULL, 
    ...)

grid.plot.SpliceExprSet(x, probes.opt = list(), expr.opt = list(col = NA, lty = 1:6),
    fig.xratio = c(2, 1), fig.yratio = c(2, 1), probepos.yscale = NULL,
    ylim = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ object of \code{\link{Probes-class}},
    \code{\link{SpliceSites-class}} or
\code{\link{SpliceExprSet-class}}}
  \item{add}{ add to an existing plot }
  \item{col}{ color(s) for the probes (recycled if necessary). }
  \item{col.typeI}{ color(s) for the type I spliced out exons }
  \item{col.typeI.window}{ background color for the type I spliced
    out `windows' }
  \item{col.typeII}{ color for the type II splicing events }
  \item{expr.opt}{list of options to plot expression values}
  \item{probepos.yscale}{ specify coordinates on the y-axis for the
probes.}
  \item{probes.opt}{options to plot the probes}
  \item{fig.xratio}{ratio for the left and right parts of the plot}
  \item{fig.yratio}{ratio for the upper and lower parts of the plot}
  \item{vp}{ a viewport (grid package stuff)}
  \item{xlim}{ range for the x-axis (see \code{plot}).}
  \item{ylim}{ range for the y-axis}
  \item{\dots}{ optional parameters }
}
\details{
  The `type I'/'type II' thing is described in the references found in
  the help files for \code{\link{plot.SpliceSites}}.
}
\value{
  These functions are mainly used for their side effects.
  \code{grid.plot.SpliceSites} returns the range for the y-axis when
  needed.
}
\seealso{ plot.SpliceSites, plot.Probes}
\examples{
## plot splice sites
data(spsites)
grid.plot(spsites)

## plot probes
data(probes)
grid.plot(probes)

## combined plot
grid.plot(probes, spsites)

}
\keyword{ hplot }
