\name{plot.SpliceSites}
\alias{plot.SpliceSites}
\alias{plot.Probes}
\title{ plot a SpliceSites object }
\description{
  plot objects.
}
\usage{
plot.Probes(x, col="black",
            xlab = "sequence", ylab = "probes",   
            add=FALSE, probepos.yscale=NULL, xlim=NULL,
            ...)
plot.SpliceSites(x, col.typeI = "orange",
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
  \item{xlim, ylim}{range of plotting window}
  \item{xlab, ylab}{labels for the axis}
  \item{\dots}{ optional parameters to be passed to the function \code{plot}. }
}
\details{
  If the parameter \code{main} is not specified, the function tries to
  extract the attribute `name' from \code{x}. 
  
  The two functions can be combined to display both objects on the same plot.
}
\value{
  The range for the y-axis is returned whenever needed (see \code{invisible}).
}
\references{"Standardized output for putative alternative splicing; a R package as an 
application to combine them with microarray data", Gautier L. Dao C. and
Yang U.C., 2003, submitted.}
\author{ Laurent }

\seealso{\code{\link{SpliceSites-class}}}
\examples{
data(spsites)

plot(spsites, main=attr(spsites, "name"))

sp.pData <- spsites@spsiteIpos.pData

##col <- as.integer(factor(sp.pData$tissue))

##plot(spsites, col.typeI=col, main=attr("name", spsites))

}
\keyword{ hplot }% at least one, from doc/KEYWORDS
