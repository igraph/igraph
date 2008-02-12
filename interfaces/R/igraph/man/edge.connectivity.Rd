\name{edge.connectivity}
\alias{edge.connectivity}
\alias{edge.disjoint.paths}
\alias{graph.adhesion}
\concept{Edge connectivity}
\concept{Edge-disjoint paths}
\concept{Graph adhesion}
\title{Edge connectivity.}
\description{The edge connectivity of a graph or two vertices, this is
  recently also called group adhesion.}
\usage{
edge.connectivity(graph, source=NULL, target=NULL, checks=TRUE)
edge.disjoint.paths(graph, source, target)
graph.adhesion(graph, checks=TRUE)
}
\arguments{
  \item{graph}{The input graph.}
  \item{source}{The id of the source vertex, for
    \code{edge.connectivity} it can be \code{NULL}, see details below.}
  \item{target}{The id of the target vertex, for
    \code{edge.connectivity} it can be \code{NULL}, see details below.}
  \item{checks}{Logical constant. Whether to check that the graph is
    connected and also the degree of the vertices. If the graph is
    not (strongly) connected then the connectivity is obviously zero. Otherwise
    if the minimum degree is one then the edge connectivity is also
    one. It is a good idea to perform these checks, as they can be
    done quickly compared to the connectivity calculation itself. 
    They were suggested by Peter McMahan, thanks Peter.
  }
}
\details{
  The edge connectivity of a pair of vertices (\code{source} and
  \code{target}) is the minimum number of edges needed to remove to
  eliminate all (directed) paths from \code{source} to \code{target}.
  \code{edge.connectivity} calculates this quantity if both the
  \code{source} and \code{target} arguments are given (and not
  \code{NULL}). 

  The edge connectivity of a graph is the minimum of the edge
  connectivity of every (ordered) pair of vertices in the graph.
  \code{edge.connectivity} calculates this quantity if neither the
  \code{source} nor the \code{target} arguments are given (ie. they are
  both \code{NULL}).

  A set of edge disjoint paths between two vertices is a set of paths
  between them containing no common edges. The maximum number of edge
  disjoint paths between two vertices is the same as their edge
  connectivity.

  The adhesion of a graph is the minimum number of edges needed to
  remove to obtain a graph which is not strongly connected. This is the
  same as the edge connectivity of the graph.

  The three functions documented on this page calculate similar
  properties, more precisely the most general is
  \code{edge.connectivity}, the others are included only for having more
  descriptive function names.
}
\value{
  A scalar real value.
}
\references{Douglas R. White and Frank Harary: The cohesiveness of
  blocks in social networks: node connectivity and conditional
  density, TODO: citation}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{graph.maxflow}}, \code{\link{vertex.connectivity}},
  \code{\link{vertex.disjoint.paths}}, \code{\link{graph.cohesion}}}
\examples{
g <- barabasi.game(100, m=1)
g2 <- barabasi.game(100, m=5)
edge.connectivity(g, 99, 0)
edge.connectivity(g2, 99, 0)
edge.disjoint.paths(g2, 99, 0)

g <- erdos.renyi.game(50, 5/50)
g <- as.directed(g)
g <- subgraph(g, subcomponent(g, 1))
graph.adhesion(g)
}
\keyword{graphs}
