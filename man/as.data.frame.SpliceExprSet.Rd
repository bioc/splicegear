\name{as.data.frame.SpliceExprSet}
\alias{as.data.frame.SpliceExprSet}
\title{ convert to a data.frame }
\description{
  Convert a SpliceExprSet object to a data.frame
}
\usage{
as.data.frame.SpliceExprSet(x, row.names = NA, optional = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object \code{\link{SpliceExprSet-class}} }
  \item{row.names}{ `NULL' or a character vector giving the row names for the
    data frame.  Missing values are not allowed}
  \item{optional}{ logical. If `TRUE', setting row names is optional.}
}
\details{
  Data are traditionally stored in objects of class \code{data.frame}.
  This function links the object-oriented design of the package with
  the large amount of functions working on data.frames.
}
\value{
A \code{data.frame}. The first variable names are \code{begin},
\code{end}, \code{isintypeI}, \code{isintypeII}, \code{exprs} and
\code{genenames}. The next variable names will the ones in the
\code{\link[Biobase]{phenoData-class} attribute of the
  \code{\link[Biobase]{exprSet-class}} object belonging the
  \code{\link{SpliceExprSet-class}}. The last variable names will be the
  ones in the slot \code{info} of the \code{\link{Probes-class}} object.}
}
\author{ Laurent }

\seealso{ }
\examples{
data(spliceset.Hs.44)

dataf <- as.data.frame(spliceset.Hs.44)

lm.panel <- function(x, y, ...) {
                                  points(x,y,...)
                                  p.lm <- lm(y~x); abline(p.lm)
                                }

## probe intensity values conditioned by the position of the probes on
## the mRNA
coplot(log(exprs) ~ Material | begin, data=dataf, panel=lm.panel)

}
\keyword{ manip }

