\name{graph.structure}
\alias{add.edges}
\alias{add.vertices}
\alias{delete.edges}
\alias{delete.vertices}
\title{Method for structural manipulation of graphs}
\description{These are the methods for simple manipulation of graphs:
  adding and deleting edges and vertices.}
\usage{
add.edges(graph, edges)
add.vertices(graph, nv)
delete.edges(graph, edges)
delete.vertices(graph, v)
}
\arguments{
  \item{graph}{The graph to work on.}
  \item{edges}{Numeric vector, the edges to add or remove.}
  \item{nv}{Numeric constant, the number of vertices to add.}
  \item{v}{Numeric vector, the vertices to remove.}
}
\details{
  \code{add.edges} adds the specified edges to the graph. The ids of the
  vertices are preserved.

  \code{add.vertices} adds the specified number of isolated vertices to
  the graph. The ids of the old vertices are preserved.

  \code{delete.edges} removes the specified edges from the graph. If a
  specified edge is not present, the function gives an error message.
  The ids of the vertices are preserved.

  \code{delete.vertices} removes the specified vertices from the graph
  together with their adjacent edges. The ids of the vertices are
  \emph{not} present.
}
\value{The new graph.}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
g <- graph.ring(10)
add.edges(g, c(1,5,2,6) )
delete.edges(g, c(0,9, 1,2) )
delete.vertices(g, c(1,6,7) )
}
\keyword{graphs}
