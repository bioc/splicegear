\name{SpliceExprSet-class}
\docType{class}
\alias{SpliceExprSet-class}
\alias{spliceSites}
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
    \item{exprs}{\code{signature(object = "SpliceExprSet")}: an accessor
    to access directly the slot \code{exprs} of the attribute \code{eset}}
    \item{plot}{\code{signature(x = "SpliceExprSet", y = "missing")}: a
      ploting method.}
    \item{show}{\code{signature(object = "SpliceExprSet")}: a printing method. }
    \item{spliceSites}{\code{signature(object = "SpliceExprSet")}:
      accessor. }
  }
}
\references{ a manuscript in preparation } %FIXME: poster
\author{ laurent@cbs.dtu.dk}

\seealso{ \code{\link{as.data.frame.SpliceExprSet}}, \code{\link{sort.SpliceExprSet}}
 and \code{\link{SpliceSites-class}} }

\examples{
data(eset, package="splicegear")
data(probes, package="splicegear")
data(spsites, package="splicegear")

spliceset <- new("SpliceExprSet", eset=eset,
                       probes=probes, spliceSites=spsites)

plot(spliceset)
}
\keyword{classes}
