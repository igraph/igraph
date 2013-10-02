\name{graph.complementer}
\alias{graph.complementer}
\concept{Graph operators}
\title{Complementer of a graph}
\description{A complementer graph contains all edges that were not
  present in the input graph.}
\usage{
graph.complementer(graph, loops=FALSE)
}
\arguments{
  \item{graph}{The input graph, can be directed or undirected.}
  \item{loops}{Logical constant, whether to generate loop edges.}
}
\details{
  \code{graph.complementer} creates the complementer of a graph. Only
  edges which are \emph{not} present in the original graph will be
  included in the new graph.

  \code{graph.complementer} keeps graph and vertex attriubutes, edge
  attributes are lost.  
}
\value{
  A new graph object.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% \seealso{}
\examples{
## Complementer of a ring
g <- graph.ring(10)
graph.complementer(g)

## A graph and its complementer give together the full graph
g <- graph.ring(10)
gc <- graph.complementer(g)
gu <- graph.union(g, gc)
gu
graph.isomorphic(gu, graph.full(vcount(g)))
}
\keyword{graphs}
