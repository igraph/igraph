\name{get.edge.ids}
\alias{get.edge.ids}
\title{Find the edge ids based on the incident vertices of the edges}
\description{
  Find the edges in an igraph graph that have the specified end
  points. This function handles multi-graph (graphs with multiple edges)
  and can consider or ignore the edge directions in directed graphs.
}
\usage{
get.edge.ids(graph, vp, directed = TRUE, error = FALSE, multi = FALSE)
}
\arguments{
  \item{graph}{The input graph.}
  \item{vp}{The indicent vertices, given as vertex ids or symbolic
    vertex names. They are interpreted pairwise, i.e. the first and
    second are used for the first edge, the third and fourth for the
    second, etc.}
  \item{directed}{Logical scalar, whether to consider edge directions in
    directed graphs. This argument is ignored for undirected graphs.}
  \item{error}{Logical scalar, whether to report an error if an edge is
    not found in the graph. If \code{FALSE}, then no error is reported,
    and zero is returned for the non-existant edge(s).}
  \item{multi}{Logical scalar, whether to handle multiple edges
    properly. If \code{FALSE}, and a pair of vertices are given twice
    (or more), then always the same edge id is reported back for
    them. If \code{TRUE}, then the edge ids of multiple edges are
    correctly reported.}
}
\details{
  igraph vertex ids are natural numbers, starting from one, up to the
  number of vertices in the graph. Similarly, edges are also numbered
  from one, up to the number of edges.

  This function allows finding the edges of the graph, via their
  incident vertices.
}
\value{A numeric vector of edge ids, one for each pair of input
  vertices. If there is no edge in the input graph for a given pair of
  vertices, then zero is reported. (If the \code{error} argument is
  \code{FALSE}.) }
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\examples{
g <- graph.ring(10)
ei <- get.edge.ids(g, c(1,2, 4,5))
E(g)[ei]

## non-existant edge
get.edge.ids(g, c(2,1, 1,4, 5,4))
}
\keyword{graphs}
