\name{getRelSignStrength}
\alias{getRelSignStrength}
\alias{getFinalRatio}
\title{ functions to perform SPLICE }
\description{
  Implementations of the SPLICE algorithm
}
\usage{
getRelSignStrength(x, tissue = as.factor(1:ncol(x)), fun = mean, nipt = 30, nitt = 30, ...)

getFinalRatio(x, tissue=as.factor(1:ncol(x)), fun=mean, ...)
}
\arguments{
  \item{x}{ a matrix. One probe per line, one column per
    sample. Typically this would be the slot \code{exprs} of an instance
    of class \code{ExprSet}.}
  \item{tissue}{ a covariate (factor) about the samples.}
  \item{fun}{ a function to obtain a summary value (\code{mean} by default) }
  \item{nipt}{ see reference. }
  \item{nitt}{ see reference. }
  \item{\dots}{ optional parameters for the function \code{fun} }
}
\details{
  \code{getFinalRatio} will call \code{getRelSignStrength}. The
  values are log-transformed. It is probably a good idea to avoid
  feeding function with values that are already on log scale.
}
\value{
  A matrix of the same dimension than the input \code{x}, holding
  'RSS' (Relative Signal Strength) or 'final ratios' respectively, as described in the reference. Two
  attributes \code{nip} and \code{nit} are attached the returned matrix.
}
\references{ Genome Research (2001), Hu et. al., vol. 11, p.1244 }
\author{ laurent@cbs.dtu.dk }
\examples{
data(spliceset)

## The intensity values in the example are log-transformed.
## Undo by taking the exponential
exprs(spliceset) <- exp(exprs(spliceset))

## Re-order the rows of different slots to have the probes sorted by
## position
spliceset <- sort.SpliceExprSet(spliceset)
## extract the expression matrix
expr.m <- exprs(spliceset)
fr <- getFinalRatio(expr.m, tissue=pData(spliceset@eset)[[1]])

}
\keyword{ manip }
