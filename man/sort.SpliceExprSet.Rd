\name{sort.SpliceExprSet}
\alias{sort.SpliceExprSet}
\title{ A function to sort a SpliceExprSet }
\description{
  Sort the probes in a SpliceExprSet (and reflect this in all the
  relevant places.
}
\usage{
sort.SpliceExprSet(x, decreasing, fun = function(x) order(x@probes@pos[, 1]), reverse = FALSE, ...)
}
\arguments{
  \item{x}{ a \code{SpliceExprSet}.}
  \item{decreasing}{ currently ignored }
  \item{fun}{ a function to do the sorting }
  \item{reverse}{ return the reverse of the sorting order }
  \item{\dots}{currently ignored}
}
\details{

}
\value{
  An object of class \code{SpliceExprSet}
}
\author{ Laurent }
\seealso{\code{\link{SpliceExprSet-class}}}
\examples{
data(spliceset)

s.spliceset <- sort.SpliceExprSet(spliceset)

}
\keyword{ manip }
