\name{grid.plot.Probes}
\alias{grid.plot.Probes}
\alias{grid.plot.SpliceSites}
%- Also NEED an '\alias' for EACH other topic documented here.
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

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ object of \code{\link{Probes-class}} or
    \code{\link{SpliceSites-class}} }
  \item{col}{ color(s) for the probes (recycled if necessary). }
  \item{col.typeI}{ color(s) for the type I spliced out exons }
  \item{col.typeI.window}{ background color for the type I spliced
    out `windows' }
  \item{col.typeII}{ color for the type II splicing events }
  \item{add}{ add to an existing plot }
  \item{probepos.yscale}{ specify coordinates on the y-axis for the probes.}
  \item{xlim}{ range for the x-axis (see \code{plot}).}
  \item{ylim}{ range for the y-axis}
  \item{vp}{ a viewport (grid package stuff)}
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
