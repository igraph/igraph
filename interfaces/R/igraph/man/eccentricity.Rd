\name{eccentricity}
\alias{eccentricity}
\alias{radius}
\concept{Eccentricity}
\concept{Radius}
\title{Eccentricity and radius}
\description{The eccentricity of a vertex is its shortest path distance
  from the farthest other node in the graph. The smallest eccentricity
  in a graph is called its radius}
\usage{
eccentricity(graph, vids=V(graph), mode=c("all", "out", "in", "total"))
radius(graph, mode=c("all", "out", "in", "total"))
}
\arguments{
  \item{graph}{The input graph, it can be directed or undirected.}
  \item{vids}{The vertices for which the eccentricity is calculated.}
  \item{mode}{Character constant, gives whether the shortest paths to or
    from the given vertices should be calculated for directed graphs. If
    \code{out} then the shortest paths \emph{from} the vertex, if
    \code{in} then \emph{to} it will be considered. If \code{all}, the
    default, then the corresponding undirected graph will be used, edge
    directions will be ignored. This argument is ignored for
    undirected graphs.}
}
\details{
 The eccentricity of a vertex is calculated by measuring the shortest
 distance from (or to) the vertex, to (or from) all vertices in the
 graph, and taking the maximum.

 This implementation ignores vertex pairs that are in different
 components. Isolate vertices have eccentricity zero.  
}
\value{
  \code{eccentricity} returns a numeric vector, containing the
  eccentricity score of each given vertex.

  \code{radius} returns a numeric scalar.
}
\references{
  Harary, F. Graph Theory. Reading, MA: Addison-Wesley, p. 35, 1994.
}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{shortest.paths}} for general shortest path
  calculations.}
\examples{
g <- graph.star(10, mode="undirected")
eccentricity(g)
radius(g)
}
\keyword{graphs}
