\name{plot.SpliceSites}
\alias{plot.SpliceSites}
\alias{plot.Probes}
\title{ plot a SpliceSites object }
\description{
  plot objects.
}
\usage{
plot.Probes(x, col="black", add=FALSE,
            probepos.yscale=NULL, xlim=NULL, ...)
plot.SpliceSites(x, col.typeI = "greenyellow",
            col.typeI.window = "yellow",
            col.typeII = "red",
            add=FALSE, ylim=NULL, ...)
}
\arguments{
  \item{x}{ object of class \code{Probes} or \code{SpliceSites}. }
  \item{col}{ color argument for the probes. }
  \item{col.typeI}{ color argument for the type I splice sites }
  \item{col.typeI.window}{ color argument for the type I `window' }
  \item{col.typeII}{ color argument for the type II splice sites }
  \item{add}{ add the plot to an existing plot. Make a new plot if `FALSE'}
  \item{probepos.yscale}{ scaling argument }
  \item{fig.yratio}{ a vector of length 2 giving the ratio of the 2
    parts of the plot (see details). }
  \item{xlim}{range of plotting window}
  \item{ylim}{range of plotting window}
  \item{\dots}{ optional parameters to be passed to the function \code{plot}. }
}
\details{
  If the parameter \code{main} is not specified, the function tries to
  extract the attribute `name' from \code{x}. 
  
  The two functions can be combined to display both objects on the same plot.
}
\value{
  The range of $y$ values.
}
\references{"Standardazed output for putative alternative splicing: application to the study of 
splice variants with microarrays", Gautier L. et al., 2003, manuscript
in preparation.}
\author{ Laurent }

\seealso{\code{\link{SpliceSites-class}}}
\examples{
data(spsites.Hs.44)

plot(spsites.Hs.44, main="Hs.44")

sp.pData <- spsites.Hs.44@spsiteIpos.pData

##col <- as.integer(factor(sp.pData$tissue))

##plot(spsites.Hs.44, col.typeI=col, main="Hs.44")

}
\keyword{ hplot }% at least one, from doc/KEYWORDS
