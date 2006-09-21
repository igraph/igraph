\name{graph.structure}
\alias{add.edges}
\alias{add.vertices}
\alias{delete.edges}
\alias{delete.vertices}
\title{Method for structural manipulation of graphs}
\description{These are the methods for simple manipulation of graphs:
  adding and deleting edges and vertices.}
\usage{
add.edges(graph, edges, \dots, attr=list())
add.vertices(graph, nv, \dots, attr=list())
delete.edges(graph, edges)
delete.vertices(graph, v)
}
\arguments{
  \item{graph}{The graph to work on.}
  \item{edges}{Edge sequence, the edges to remove.}
  \item{\dots}{Additional parameters will be added as edge/vertex
    attributes. Note that these arguments have to be named.}
  \item{attr}{Additional edge/vertex attributes to add. This will be
    concatenated to the other supplied attributes.}
  \item{nv}{Numeric constant, the number of vertices to add.}
  \item{v}{Vector sequence, the vertices to remove.}
}
\details{
  \code{add.edges} adds the specified edges to the graph. The ids of the
  vertices are preserved. The additionally supplied named arguments will
  be added as edge attributes for the new edges. If an attribute was not
  present in the original graph, its value for the original edges will
  be \code{NA}.

  \code{add.vertices} adds the specified number of isolate vertices to
  the graph. The ids of the old vertices are preserved. The additionally
  supplied named arguments will be added as vertex attributes for the
  new vertices. If an attribute was not present in the original graph,
  its value is set to \code{NA} for the original vertices.

  \code{delete.edges} removes the specified edges from the graph. If a
  specified edge is not present, the function gives an error message,
  and the original graph remains unchanged.
  The ids of the vertices are preserved.

  \code{delete.vertices} removes the specified vertices from the graph
  together with their adjacent edges. The ids of the vertices are
  \emph{not} preserved.
}
\value{The new graph.}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
g <- graph.ring(10)
add.edges(g, c(1,5,2,6) )
delete.edges(g, E(g, P=c(0,9, 1,2)) )
delete.vertices(g, c(1,6,7) )
}
\keyword{graphs}
