\name{matchprobes2Probes}
\alias{matchprobes2Probes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ create Probes object from matchprobes results }
\description{
  Create Probes object from results the results of the function
  matchprobes (in the package 'matchprobes').
}
\usage{
matchprobes2Probes(mpo, probes.length, names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mpo}{\code{Probes-class} object}
  \item{probes.length}{ Length for the probes (see details). }
  \item{names}{ names for the elements in the list returned.} 
}
\details{
  Currently only probes of unique length are assumed.
  In the case of Affymetrix chips, 25 base pairs is the value you
  probably want.
}
\value{
  A \code{list} of \code{Probes-class} objects.
}
\references{ http://www.cbs.dtu.dk/laurent/download/splicegear/ }

\seealso{ \code{\link{Probes-class}}, the package \code{matchprobes}}

\keyword{ hplot }
