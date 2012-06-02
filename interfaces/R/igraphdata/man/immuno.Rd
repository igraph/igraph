\name{immuno}
\docType{data}
\alias{immuno}
\title{Immunoglobulin interaction network}
\description{
  The undirected and connected network of interactions
  in the immunoglobulin protein. It is made up of 1316 vertices
  representing amino-acids and an edge is drawn between two
  amino-acids if the shortest distance between their C_alpha atoms
  is smaller than the threshold value \eqn{\theta=8}{theta=8} Angstrom.
}
\usage{
immuno
}
\format{
  An undirected \code{igraph} graph object.

  Graph attributes: \sQuote{name}, \sQuote{Citation}, \sQuote{Author}.
}
\source{
  See reference below.
}
\references{
  D. Gfeller, Simplifying complex networks: from a clustering to a
  coarse graining strategy, \emph{PhD Thesis EPFL}, no 3888, 2007.
  \url{http://library.epfl.ch/theses/?nr=3888}
}
\keyword{datasets}
