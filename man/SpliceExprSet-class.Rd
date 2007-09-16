\name{SpliceExprSet-class}
\docType{class}
\alias{SpliceExprSet-class}

\alias{grid.plot}
\alias{spliceSites}
%\alias{grid.plot.SpliceExprSet}
\alias{grid.plot,SpliceExprSet,missing-method}
\alias{plot,SpliceExprSet,missing-method}
\alias{show,SpliceExprSet-method}
\alias{spliceSites,SpliceExprSet-method}
\alias{exprs,SpliceExprSet-method}
\alias{exprs<-,SpliceExprSet,ANY-method}
\alias{grid.plot,SpliceSites,missing-method}


\title{Class "SpliceExprSet"}
\description{A class to store probe expression values with alternative
  splicing information}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SpliceExprSet", ...)}.

}
\section{Slots}{
  \describe{
    \item{\code{spliceSites}:}{Object of class \code{"SpliceSites"}. The
    probes and splice site information.}
    \item{\code{probes}:}{Object of class \code{"Probes"}. The matching expression
      values.}
    \item{\code{eset}:}{Object of class \code{"ExpressionSet"}. The matching expression
      values.}
  }
}
\section{Methods}{
  \describe{
    \item{grid.plot}{\code{signature(x = "SpliceExprSet", y = "missing")}: ... }
    \item{plot}{\code{signature(x = "SpliceExprSet", y = "missing")}: a
      ploting method.}
    \item{show}{\code{signature(object = "SpliceExprSet")}: a printing method. }
    \item{spliceSites}{\code{signature(object = "SpliceExprSet")}:
      accessor. }
  }
}
\references{ a manuscript in preparation } %FIXME: poster
\author{ laurent@cbs.dtu.dk}

\seealso{\code{\link{as.data.frame.SpliceExprSet}}, \code{\link{sort.SpliceExprSet}}
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
