\name{SpliceExprSet-class}
\docType{class}
\alias{SpliceExprSet-class}
\title{Class "SpliceExprSet" ~~~ }
\description{ A class to store probe expression values with alternative
  splicing information}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SpliceExprSet", ...)}.
    ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{spliceSites}:}{Object of class \code{"SpliceSites"}. The
    probes and splice site information.}
    \item{\code{eset}:}{Object of class \code{"exprSet"}. The matching expression
    values.}
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "SpliceExprSet")}: a printing method. }
    \item{exprs}{\code{signature(object = "SpliceExprSet")}: an accessor
    to access directly the slot \code{exprs} of the attribute \code{eset}}
    \item{plot}{\code{signature(x = "SpliceExprSet", y = "missing")}: a
      ploting method.}
    \item{spliceSites}{\code{signature(object = "SpliceExprSet")}: accessor.}
  }
}
\references{ a manuscript in preparation } %FIXME: poster
\author{ laurent@cbs.dtu.dk}

\seealso{ \code{\link{as.data.frame.SpliceExprSet}}, \code{\link{sort.SpliceExprSet}}
 and \code{\link{SpliceSites-class}} }

\examples{
data(eset.Hs.44)
data(probes.Hs.44)
data(spsites.Hs.44)

spliceset.Hs.44 <- new("SpliceExprSet", eset=eset.Hs.44,
                       probes=probes.Hs.44, spliceSites=spsites.Hs.44)

plot(spliceset.Hs.44)
}
\keyword{classes}
