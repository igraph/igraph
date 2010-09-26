\name{biconnected.components}
\alias{biconnected.components}
\concept{Biconnected component}
\title{Biconnected components}
\description{Finding the biconnected components of a graph}
\usage{
biconnected.components(graph)
}
\arguments{
  \item{graph}{The input graph. It is treated as an undirected graph,
    even if it is directed.}
}
\details{
  A graph is biconnected if the removal of any single vertex (and its
  adjacent edges) does not disconnect it.

  A biconnected component of a graph is a maximal biconnected subgraph
  of it. The biconnected components of a graph can be given by the
  partition of its edges: every edge is a member of exactly one
  biconnected component. Note that this is not true for vertices: the same
  vertex can be part of many biconnected components.  
}
\value{
  A named list with three components:
  \item{no}{Numeric scalar, an integer giving the number of
    biconnected components in the graph.}
  \item{tree_edges}{The components themselves, a list of numeric
    vectors. Each vector is a set of edge ids giving the edges in a
    biconnected component. These edges define a spanning tree of the
    component.}
  \item{component_edges}{A list of numeric vectors. It gives all edges
  in the components.}
  \item{components}{A list of numeric vectors, the vertices of the
    components.}
  \item{articulation_points}{The articulation points of the graph. See
    \code{\link{articulation.points}}.}
}
%\references
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{articulation.points}}, \code{\link{clusters}},
  \code{\link{is.connected}}, \code{\link{vertex.connectivity}}}
\examples{
g <- graph.disjoint.union( graph.full(5), graph.full(5) )
clu <- clusters(g)$membership
g <- add.edges(g, c(which(clu==1), which(clu==2)))
bc <- biconnected.components(g)
}
\keyword{graphs}
