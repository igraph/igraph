\name{neighborhood}
\alias{neighborhood}
\alias{neighborhood.size}
\alias{graph.neighborhood}
\alias{connect.neighborhood}
\title{Neighborhood of graph vertices}
\description{These functions find the vertices not farther than a given
  limit from another fixed vertex, these are called the neighborhood of
  the vertex.}
\usage{
neighborhood.size(graph, order, nodes=V(graph), mode=c("all", "out", "in"))
neighborhood(graph, order, nodes=V(graph), mode=c("all", "out", "in"))
graph.neighborhood(graph, order, nodes=V(graph), mode=c("all", "out", "in"))
connect.neighborhood(graph, order, mode=c("all", "out", "in", "total"))
}
\arguments{
  \item{graph}{The input graph.}
  \item{order}{Integer giving the order of the neighborhood.}
  \item{nodes}{The vertices for which the calculation is performed.}
  \item{mode}{Character constatnt, it specifies how to use the direction
    of the edges if a directed graph is analyzed. For \sQuote{out} only
    the outgoing edges are followed, so all vertices reachable from the
    source vertex in at most \code{order} steps are counted. For
    \sQuote{"in"} all vertices from which the source vertex is reachable
    in at most \code{order} steps are counted. \sQuote{"all"} ignores
    the direction of the edges. This argument is ignored for undirected
    graphs.} 
}
\details{
  The neighborhood of a given order \code{o} of a vertex \code{v}
  includes all vertices which are closer to \code{v} than the
  order. Ie. order 0 is always \code{v} itself, order 1 is \code{v}
  plus its immediate neighbors, order 2 is order 1 plus the immediate
  neighbors of the vertices in order 1, etc.

  \code{neighborhood.size} calculates the size of the neighborhoods for
  the given vertices with the given order.

  \code{neighborhood} calculates the neighborhoods of the given vertices
  with the given order parameter.

  \code{graph.neighborhood} is creates (sub)graphs from all neighborhoods of
  the given vertices with the given order parameter. This function
  preserves the vertex, edge and graph attributes.

  \code{connect.neighborhood} creates a new graph by connecting each
  vertex to all other vertices in its neighborhood.
}
\value{
  \code{neighborhood.size} returns with an integer vector.

  \code{neighborhood} returns with a list of integer vectors.

  \code{graph.neighborhood} returns with a list of graphs.

  \code{connect.neighborhood} returns with a new graph object.
}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}, the first version was
  done by Vincent Matossian}
\examples{
g <- graph.ring(10)
neighborhood.size(g, 0, 1:3)
neighborhood.size(g, 1, 1:3)
neighborhood.size(g, 2, 1:3)
neighborhood(g, 0, 1:3)
neighborhood(g, 1, 1:3)
neighborhood(g, 2, 1:3)

# attributes are preserved
V(g)$name <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j")
graph.neighborhood(g, 2, 1:3)

# connecting to the neighborhood
g <- graph.ring(10)
g <- connect.neighborhood(g, 2)
}
\keyword{graphs}
