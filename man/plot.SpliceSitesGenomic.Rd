\name{plot.SpliceSitesGenomic}
\alias{plot.SpliceSitesGenomic}
\title{ Function to plot SpliceSitesGenomic objects }
\description{
  Function to plot SpliceSitesGenomic objects }
}
\usage{
plot.SpliceSitesGenomic(x, col.variant = par("col"), col.exon = "white",
                        split = FALSE, main = NULL, ...)
}
\arguments{
  \item{x}{ \code{SpliceSitesGenomic-class} }
  \item{col.variant}{ a vector of colors for the different variants. The
    colors are recycled as necessary.}
  \item{col.exon} { a vector of colors for the exons. The
    colors are recycled as necessary.}
  \item{split}{ split the plot of the variants in individual plots}
  \item{main}{\code{character} to use as a title. Recycled as necessary.}
  \item{\dots}{ optional graphical parameters }
}
\value{
  This function is used for its side-effect.
}
\author{ Laurent }
\seealso{ \code{\link{SpliceSitesGenomic-class}} }
\examples{
## a 10 bp window
seq.length <- as.integer(10)
## positions of the exons
spsiteIpos <- matrix(c(1, 3.5, 5, 9, 3, 4, 8, 10), nc=2)
## known variants
variants <- list(a=c(1,2,3,4), b=c(1,2,3), c=c(1,3,4))
##
n.exons <- nrow(spsiteIpos)

spvar <- new("SpliceSitesGenomic", spsiteIpos=spsiteIpos,
         variants=variants, seq.length=seq.length)

par(mfrow = c(3,1), mar = c(3.1, 2.1, 2.1, 1.1))

plot(spvar, split=TRUE, col.exon=rainbow(n.exons))

}
\keyword{ hplot }
