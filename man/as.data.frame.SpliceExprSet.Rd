\name{as.data.frame.SpliceExprSet}
\alias{as.data.frame.SpliceExprSet}
\alias{as.data.frame.SpliceSites}
\title{SpliceExprSet object to data.frame converter}
\description{
  Converts a SpliceExprSet object to a data.frame.
}
\usage{
\method{as.data.frame}{SpliceSites}(x, row.names = NA, optional = NA, \dots)

\method{as.data.frame}{SpliceExprSet}(x, row.names = NA, optional = NA, \dots)
}
\arguments{
  \item{x}{object \code{\link{SpliceSites-class}} or \code{\link{SpliceExprSet-class}}.}
  \item{row.names}{\code{NULL} or a character vector giving the row names for the
    data frame.  Missing values are not allowed.}
  \item{optional}{logical. If \code{TRUE}, setting row names is optional.}
  \item{\dots}{currently ignored.}
}
\details{
  Data are traditionally stored in objects of class \code{data.frame}.
  This function links the object-oriented design of the package with
  the large amount of functions working on data.frames.
}
\value{
A \code{data.frame}. For both functions the first column names are \code{begin},
\code{end}, \code{isintypeI}, \code{isintypeII}, \code{exprs} and \code{genenames}.
In the case of \code{as.data.frame.SpliceExprSet}, the next variable names will the ones in the
\code{\link[Biobase]{AnnotatedDataFrame-class}} attribute of the
\code{\link[Biobase]{ExpressionSet-class}} object belonging the
\code{\link{SpliceExprSet-class}}. The last variable names will be the
ones in the slot \code{info} of the \code{\link{Probes-class}} object.
}
\author{Laurent Gautier}
\examples{
data(spliceset)

dataf <- as.data.frame(spliceset)

lm.panel <- function(x, y, ...) {
                                  points(x,y,...)
                                  p.lm <- lm(y~x); abline(p.lm)
                                }

## probe intensity values conditioned by the position of the probes on
## the mRNA
coplot(log(exprs) ~ Material | begin, data=dataf, panel=lm.panel)
}
\keyword{manip}
