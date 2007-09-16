\name{barplot.SpliceSites}
\alias{barplot.SpliceSites}
\title{ barplot for SpliceSites }
\description{
  Display a barplot of the associated phenoData
}
\usage{
\method{barplot}{SpliceSites}(height, type.as = c("typeI", "typeII", "all"),
                    info = "tissue", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{height}{ object of class \code{\link{SpliceSites-class}} }
  \item{type.as}{ the type of alternative splicing (see
    \code{SpliteSites-class} for futher details) }
  \item{info}{ the name of the covariate in the phenoData (see details)}
  \item{\dots}{ optional parameters to be passed to the underlying
    function \code{\link{barplot}}. }
}
\details{
  When the objects are built from the XML format we propose as an
  exchange, the parameter \code{info} can at least take the values
  \code{"tissue"} and \code{"histology"}. One can refer to the
  slots \code{spsiteIpos.pData} and \code{spsiteIIpos.pData} to
  know what are the possible choices.
}
\value{
  See the value returned by the function \code{\link{barplot}}.
}
\seealso{ \code{\link{SpliceSites-class}}, \code{\link{barplot}} }
\examples{

data(spsites)

barplot(spsites)

}
\keyword{ hplot }
