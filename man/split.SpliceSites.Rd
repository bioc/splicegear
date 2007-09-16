\name{split.SpliceSites}
\alias{split.SpliceSites}

\title{ split an instance of SpliceSites }
\description{
  Split an instance of SpliceSites into several instances of SpliceSites
}
\usage{
\method{split}{SpliceSites}(x, f = list(typeI = NA, typeII = NA), drop=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an instance of class \code{\link{SpliceSites-class}}.  }
  \item{f}{ a list of two factors (see details). }
  \item{drop}{not used (here to keep R CMD check happy)}
  \item{...}{see \code{drop} above.}
}
\details{
  The split usually is performed on a factor. Two factors are required
  because of the two categories of splicing events (type I (deletion of
  a fragment of the reference sequence) and type II (insertion of an
  element of the reference sequence). A \code{character} can be used
  instead of a \code{factor}. In this case the covariates with the
  given name, in the slots \code{spsiteIpos.pData} and
  \code{spsiteIIpos.pData} respectively, are used to make the split.
  When equal to \code{NA}, the covariate named \code{site} will be
  used.
}
\value{
  A list of objects of class \code{SpliceSites}
}
\seealso{ \code{\link{SpliceSites-class}} }
\examples{
data(spsites)
split(spsites)
}
\keyword{ manip }
