\name{buildSpliceSites}
\alias{buildSpliceSites}
\alias{queryPALSdb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Functions to query PALSdb }
\description{
  Functions to make a query on PALSdb, and build objects from the result
  of a query.
}
\usage{
queryPALSdb(query, disp = c("data", "browser"),
            field=c("keyword", "ug_id", "gb_id", "cluster_count"))
buildSpliceSites(xml)
}
\arguments{
  \item{query}{query string}
  \item{xml}{an object of class XML (as returned by
    \code{\link[XML]{xmlTreeParse}})}
  \item{disp}{(idem \code{\link[annotate]{genbank}} and \code{\link[annotate]{pubmed}})}
  \item{field}{The field on which the query will be based}
}
\details{
  \code{queryPALSdb} returns an an object of class XML when \code{disp =
    "data"}.
}
\value{
  
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
