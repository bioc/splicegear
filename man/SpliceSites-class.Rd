\name{SpliceSites-class}
\docType{class}
\alias{SpliceSites-class}
\title{Class "SpliceSites"}
\description{ A class to store (putative) splice sites }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SpliceSites", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{probepos}:}{Object of class \code{"matrix"}. A two-columns
      matrix to store the begin and end positions of probes (one row per
      probe).}
    \item{\code{seq}:}{Object of class \code{"character"}. The reference
      sequence.}
    \item{\code{seq.length}:}{Object of class \code{"integer"}. The length for
      the reference sequence (used when the slot \code{seq} is set to
      \code{""}).}
    \item{\code{spsiteIpos}:}{Object of class \code{"matrix"}. A
      two-columns
      matrix to store the begin and end positions of type I splice variant.}
    \item{\code{spsiteIIpos}:}{Object of class \code{"integer"}. A
      vector to store the positions for type II splice variants.}
    \item{\code{spsiteIIIpos}:}{Object of class \code{"matrix"}. Idem
      \code{spsiteIpos}, but for type III splice variants.
    }
    \item{\code{spsiteIpos.pData}:}{Object of class
      \code{\link[Biobase]{"phenoData"}}. Used to store covariate
      information related to the splice variants}.
    \item{\code{spsiteIIpos.pData}:}{Object of class
      \code{\link[Biobase]{"phenoData"}}}.
    \item{\code{spsiteIIIpos.pData}:}{Object of class \code{\link[Biobase]{"phenoData"}}}. 
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "SpliceSites")}: A printing method. }
    \item{plot}{\code{signature(x = "SpliceSites", y = "missing")}: A
      plotting method }
  }
}
\references{ "Plenty of splicin' or 'can regular Affymetrix chips be used
  to observe alternative splicing ?'", Gautier L. et al., 2003, manuscript in
  preparation (and the title might have to chang\ldots).} %FIXME: NAR paper + poster
\author{ laurent@cbs.dtu.dk }
\seealso{ \code{\link{isSpliceSiteOnProbe}},
  \code{\link{isProbeOnSpliceSite}},
  \code{\link{plot.SpliceSites}},
  \code{\link{spliceset.Hs.44}}.
}

\examples{
data(spliceset.Hs.44)

print(spliceset.Hs.44)

par(mfrow=c(1,2))

plot(spliceset.Hs.44, main="Hs.44")

## filter out supporting matches with unique positions
filter.typeI <- function(x) {unique(x[duplicated(x), , drop=FALSE])}
spliceset.Hs.44.filter <- spliceset.Hs.44
sSites <- spliceset.Hs.44.filter@spliceSites
sSites@spsiteIpos <- filter.typeI(sSites@spsiteIpos)
spliceset.Hs.44.filter@spliceSites <- sSites
## plot the resulting new object
plot(spliceset.Hs.44.filter)

}
\keyword{classes}
