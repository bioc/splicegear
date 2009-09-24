\name{grid.expand.gp}
\alias{grid.expand.gp}
\alias{grid.make.numeric2npc}

\title{ convenience functions for grid }
\description{
  Convenience function to use the package grid
}
\usage{
grid.expand.gp(n, parlist = list())
grid.make.numeric2npc(x, xlim=NULL, lower.blank=0, upper.blank=0)

}
\arguments{
  \item{n}{ number of parameters }
  \item{parlist}{ list of parameters }
  \item{x}{numeric value}
  \item{xlim}{range for Xs}
  \item{lower.blank, upper.blank}{size for margins}
}
\details{
  call the function \code{gpar} on the list of parameters.
}
\value{
  Function used for its side effect.
}
\seealso{ \code{\link[lattice:Lattice]{lattice}}
}
\keyword{ hplot }

