\name{all_simple_paths}
\alias{all_simple_paths}
\concept{Simple path}
\title{List all simple paths from one source}
\description{This function lists are simple paths from one source vertex
  to another vertex or vertices. A path is simple if the vertices it
  visits are not visited more than once.}
\usage{
all_simple_paths (graph, from, to = V(graph),
    mode = c("out", "in", "all", "total"))
}
\arguments{
  \item{graph}{The input graph.}
  \item{from}{The source vertex.}
  \item{to}{The target vertex of vertices. Defaults to all vertices.}
  \item{mode}{Character constant, gives whether the shortest paths to or
    from the given vertices should be calculated for directed graphs. If
    \code{out} then the shortest paths \emph{from} the vertex, if
    \code{in} then \emph{to} it will be considered. If \code{all}, the
    default, then the corresponding undirected graph will be used,
    ie. not directed paths are searched. This argument is ignored for
    undirected graphs.}
}
\details{
 Note that potentially there are exponentially many
 paths between two vertices of a graph, and you may
 run out of memory when using this function, if your
 graph is lattice-like.

 This function currently ignored multiple and loop edges.
}
\value{
  A list of integer vectors, each integer vector is a path from the
  source vertex to one of the target vertices. A path is given by its
  vertex ids.
}
%\references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
%\seealso{}
\examples{
g <- make_ring(10)
all_simple_paths(g, 1, 5)
all_simple_paths(g, 1, c(3,5))
}
\keyword{graphs}
