\name{isProbeOnSpliceSite}
\alias{isProbeOnSpliceSite}
\title{ Check the presence of probes on certain exons }
\description{
  Return whether the probes are located on exons involved in (putative)
  alternative splicing or not.
}
\usage{
isProbeOnSpliceSite(probes, spSites)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probes}{ object of class \code{Probes}}
  \item{spSites}{ object of class \code{spliceSites}}
}
\details{
}
\value{
  The returned value in a list of two vectors of mode \code{logical}
  of the same length:
  \item{isintypeI}{whether the probes are in a `type I' region or not.}
  \item{isintypeII}{whether the probe are in a `type II' region or not.}
}
\references{ For details about `type I' and `type II', please refer to
  Huang Y.-H and Chen Y.-T and Lai J.-J. and Yang S.-T. and
  Yang U.-C., PALSdb: Putative Alternative Splicing database, Nucleic Acids Research,
  2002, pages 186-190 }
\author{ Laurent }
\examples{

}
\keyword{ manip }

