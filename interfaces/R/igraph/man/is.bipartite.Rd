\name{is.bipartite}
\alias{is.bipartite}
\concept{Bipartite graph}
\concept{Two-mode network}
\title{Decide whether a graph is bipartite}
\description{This function decides whether the vertices of a network can
  be mapped to two vertex types in a way that no vertices of the same
  type are connected.}
\usage{
is.bipartite(graph)
}
\arguments{
  \item{graph}{The input graph.}
}
\details{
  A bipartite graph in igraph has a \sQuote{\code{type}} vertex
  attribute giving the two vertex types.

  This function simply checks whether a graph \emph{could} be
  bipartite. It tries to find a mapping that gives a possible division
  of the vertices into two classes, such that no two vertices of the
  same class are connected by an edge.

  The existence of such a mapping is equivalent of having no circuits of
  odd length in the graph. A graph with loop edges cannot bipartite.

  Note that the mapping is not necessarily unique, e.g. if the graph has
  at least two components, then the vertices in the separate components
  can be mapped independently.
}
\value{
  A named list with two elements:
  \item{res}{A logical scalar, \code{TRUE} if the can be bipartite,
    \code{FALSE} otherwise.}
  \item{type}{A possibly vertex type mapping, a logical vector. If no
    such mapping exists, then an empty vector.}
}
% e
% \references{}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu} }
% \seealso{}
\examples{
## A ring has just one loop, so it is fine
g <- graph.ring(10)
is.bipartite(g)

## A star is fine, too
g2 <- graph.star(10)
is.bipartite(g2)

## A graph containing a triangle is not fine
g3 <- graph.ring(10)
g3 <- add.edges(g3, c(0,2))
is.bipartite(g3)
}
\keyword{graphs}
