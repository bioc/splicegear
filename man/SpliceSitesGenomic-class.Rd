\name{SpliceSitesGenomic-class}
\docType{class}
\alias{SpliceSitesGenomic-class}
\alias{plot,SpliceSitesGenomic-method}
\alias{plot,SpliceSitesGenomic,missing-method}
\title{Class "SpliceSitesGenomic"}
\description{A class to store alternative splicing information on a genomic point of view.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SpliceSitesGenomic", seq, seq.length, spsiteIpos, spsiteIIpos, spsiteIIIpos, spsiteIpos.pData, spsiteIIpos.pData, spsiteIIIpos.pData, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{variants}:}{Object of class \code{"list"}. There is one
      element per splice variant. Each element in the list should be a
      vector of integers. Each integer refers to an exon. The sequence
      of integers determines the sequence of exons in the splice variant.}
    \item{\code{seq}:}{Object of class \code{"character", from class "SpliceSites"}.}
    \item{\code{seq.length}:}{Object of class \code{"integer", from class "SpliceSites"}.}
    \item{\code{spsiteIpos}:}{Object of class \code{"matrix", from class "SpliceSites"}.}
    \item{\code{spsiteIIpos}:}{Object of class \code{"integer", from
	class "SpliceSites"}. This should not have any practical use in this class.}
    \item{\code{spsiteIIIpos}:}{Object of class \code{"matrix", from
	class "SpliceSites"}. This should not have any practical use in this class.}
    \item{\code{spsiteIpos.pData}:}{Object of class \code{"AnnotatedDataFrame", from class "SpliceSites"}.}
    \item{\code{spsiteIIpos.pData}:}{Object of class \code{"AnnotatedDataFrame",
	from class "SpliceSites"}. This should not have any practical use in this class.}
    \item{\code{spsiteIIIpos.pData}:}{Object of class \code{"AnnotatedDataFrame",
	from class "SpliceSites"}. This should not have any practical use in this class.}
  }
}
\section{Extends}{
Class \code{"SpliceSites"}, directly.
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "SpliceSitesGenomic", y =
	"missing")}: a plotting method for demonstration purposes. }
  }
}
\seealso{ \code{\link{SpliceSites-class}} and \code{\link{plot.SpliceSitesGenomic}}.}
\examples{
## a 10 bp window
seq.length <- as.integer(10)
## positions of the exons
spsiteIpos <- matrix(c(1, 3.5, 5, 9, 3, 4, 8, 10), nc=2)
## known variants
variants <- list(a=c(1,2,3,4), b=c(1,2,3), c=c(1,3,4))

spvar <- new("SpliceSitesGenomic", spsiteIpos=spsiteIpos,
         variants=variants, seq.length=seq.length)

plot(spvar)
}
\keyword{classes}
