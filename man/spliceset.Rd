\name{spliceset}
\alias{spliceset}
\alias{probes}
\alias{spsites}
\alias{eset}
\docType{data}
\title{Example data for splicegear}
\description{
  The putative splice variants for a reference sequence,
  the matching probes from the Affymetrix chip 'HG-U95A'
  and probe intensities from the 'dilution' dataset.
}
\usage{
data(eset, package="splicegear")
data(probes)
data(spsites, package="splicegear")
data(spliceset, package="splicegear")
}
\format{
  The formats are objects of class \code{\link[Biobase]{exprSet-class}},
  \code{\link{Probes-class}}, \code{\link{SpliceSites-class}} and
  \code{\link{SpliceExprSet-class}} respectively.
}
\details{
  The attribute "name" is set to the ID of the Unigene cluster from
  which the reference sequence is taken.
}
\source{
  
}
\references{
  "PALSdb", ref.
  GeneLogic's dilution dataset.
}
\examples{
data(spliceset, package="splicegear")

plot(spliceset, main=attr(spliceset, "name"))
}
\keyword{datasets}
