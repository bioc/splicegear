\name{Probes-class}
\docType{class}
\alias{Probes-class}
\alias{initialize}
\alias{show}
\alias{plot}
\alias{plot}

\alias{initialize,Probes-method}
\alias{show,Probes-method}
\alias{plot,Probes-method}
\alias{plot,Probes-method}
\alias{grid.plot,Probes,missing-method}
\alias{grid.plot,Probes,SpliceSites-method}
\alias{plot,Probes,missing-method}
\alias{plot,Probes,SpliceSites-method}

\title{Class "Probes" }
\description{Information about a set of probes }
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("Probes", pos)}
  or \code{new("Probes", pos, info)}.
  The object are primarily storing the location of the probe on a
  matching sequence sequence. Optional information can be stored
  in the slot \code{info} (a \code{data.frame}).
}
\section{Slots}{
  \describe{
    \item{\code{pos}:}{Object of class \code{"matrix"}. It expects one
      row per probe. The first column should give the start position
      while the second column should give the end position}
    \item{\code{info}:}{Object of class \code{"data.frame"}. Optional
      information one wishes to carry around can be stored here.}
  }
}

\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "Probes")}: ... }
    \item{show}{\code{signature(object)}: show minimal information}
    \item{plot}{\code{signature(x = "Probes", y = "missing")}: plot the
      position of the probes. (see \code{\link{plot.Probes}}) }
    \item{plot}{\code{signature(x = "Probes", y = "SpliceSites")}: plot
      the positions of the probes and the positions of the splice sites}
  }
}

\seealso{ \code{\link{SpliceSites-class}}, \code{\link{SpliceExprSet-class}}, }

\examples{
data(probes)

plot(probes)
}
\keyword{classes}
