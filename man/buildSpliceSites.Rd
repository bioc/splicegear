\name{buildSpliceSites}
\alias{buildSpliceSites}
\alias{queryPALSdb}
\alias{getPALSdbURL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Functions to query PALSdb }
\description{
  Functions to make a query on PALSdb, and build objects from the result
  of a query.
}
\usage{
queryPALSdb(query, disp = c("data", "browser"),
            field = c("keyword", "ug_id", "gb_id", "cluster_count"),
            species = c("human", "mouse"),
            e.value = "1e-1",
            ident.threshold = c("90% 50b", "95% 50b", "90% 45b"),
            verbose = FALSE)

getPALSdbURL(query, disp = c("data", "browser"),
            field = c("keyword", "ug_id", "gb_id", "cluster_count"),
            species = c("human", "mouse"),
            e.value = "1e-1",
            ident.threshold = c("90% 50b", "95% 50b", "90% 45b"),
            verbose = FALSE)

buildSpliceSites(xml, verbose=TRUE)
}
\arguments{
  \item{query}{query string}
  \item{xml}{an object of class XML (as returned by
    \code{\link[XML]{xmlTreeParse}})}
  \item{disp}{(idem \code{\link[annotate]{genbank}} and \code{\link[annotate]{pubmed}})}
  \item{field}{The field on which the query will be based}
  \item{species}{the specie to work with}
  \item{e.value}{}
  \item{ident.threshold}{threshold for matching sequences}
  \item{verbose}{verbose output.}
}
\details{
  \code{queryPALSdb} returns an an object of class XML when \code{disp =
    "data"}.

  The function \code{buildSpliceSites} constructs \code{SpliceSites}
  objects from the XML data. The variables in the slots
  \code{spsiteIpos.pData} and \code{spsiteIIpos.pData} are at least
  \code{tissue} (tissue information), \code{histology} and \code{site}
  (site numbering).
}
\value{
  An object of class \code{XML} for \code{queryPALSdb}, an URL for
  \code{getPALSdbURL} or a \code{list} of objects of class \code{SpliceSites}.
}
\references{"Standardazed output for putative alternative splicing: application to the study of 
splice variants with microarrays", Gautier L. et al., 2003, manuscript
in preparation.}
\author{ laurent@cbs.dtu.dk }

\seealso{\code{\link{queryPALSdb}}}
\examples{
library(XML)

filename <- system.file("data", "example.xml", package="splicegear")

xml <- xmlTreeParse(filename, asTree=TRUE)

spsites <- buildSpliceSites(xml)

}
\keyword{ manip }
