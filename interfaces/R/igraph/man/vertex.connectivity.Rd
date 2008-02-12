\name{vertex.connectivity}
\alias{vertex.connectivity}
\alias{vertex.disjoint.paths}
\alias{graph.cohesion}
\concept{Vertex connectivity}
\concept{Vertex-disjoint paths}
\concept{Graph cohesion}
\title{Vertex connectivity.}
\description{The vertex connectivity of a graph or two vertices, this is
  recently also called group cohesion.}
\usage{
vertex.connectivity(graph, source=NULL, target=NULL, checks=TRUE)
vertex.disjoint.paths(graph, source, target)
graph.cohesion(graph, checks=TRUE)
}
\arguments{
  \item{graph}{The input graph.}
  \item{source}{The id of the source vertex, for
    \code{vertex.connectivity} it can be \code{NULL}, see details
    below.}
  \item{target}{The id of the target vertex, for
    \code{vertex.connectivity} it can be \code{NULL}, see details
    below.}
  \item{checks}{Logical constant. Whether to check that the graph is
    connected and also the degree of the vertices. If the graph is
    not (strongly) connected then the connectivity is obviously zero. Otherwise
    if the minimum degree is one then the vertex connectivity is also
    one. It is a good idea to perform these checks, as they can be
    done quickly compared to the connectivity calculation itself. 
    They were suggested by Peter McMahan, thanks Peter.}
}
\details{
  The vertex connectivity of two vertices (\code{source} and
  \code{target}) in a directed graph is the minimum number of vertices
  needed to remove from the graph to eliminate all (directed) paths from
  \code{source} to \code{target}. \code{vertex.connectivity} 
  calculates this quantity if both the \code{source} and \code{target}
  arguments are given and they're not \code{NULL}.

  The vertex connectivity of a graph is the minimum vertex connectivity
  of all (ordered) pairs of vertices in the graph. In other words this
  is the minimum number of vertices needed to remove to make the graph
  not strongly connected. (If the graph is not strongly connected then
  this is zero.) \code{vertex.connectivity} calculates this quantitty if
  neither the \code{source} nor \code{target} arguments are
  given. (Ie. they are both \code{NULL}.)

  A set of vertex disjoint directed paths from \code{source} to \code{vertex}
  is a set of directed paths between them whose vertices do not contain common
  vertices (apart from \code{source} and \code{target}). The maximum number of
  vertex disjoint paths between two vertices is the same as their vertex
  connectivity. 

  The cohesion of a graph (as defined by White and Harary, see
  references), is the vertex connectivity of the graph. This is
  calculated by \code{graph.cohesion}.

  These three functions essentially calculate the same measure(s), more
  precisely \code{vertex.connectivity} is the most general, the other
  two are included only for the ease of using more descriptive function
  names.
}
\value{
  A scalar real value.
}
\references{Douglas R. White and Frank Harary: The cohesiveness of
  blocks in social networks: node connectivity and conditional
  density, TODO: citation}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{graph.maxflow}}, \code{\link{edge.connectivity}},
  \code{\link{edge.disjoint.paths}}, \code{\link{graph.adhesion}}}
\examples{
g <- barabasi.game(100, m=1)
g <- delete.edges(g, E(g)[ 99 \%--\% 0 ])
g2 <- barabasi.game(100, m=5)
g2 <- delete.edges(g2, E(g2)[ 99 \%--\% 0])
vertex.connectivity(g, 99, 0)
vertex.connectivity(g2, 99, 0)
vertex.disjoint.paths(g2, 99, 0)

g <- erdos.renyi.game(50, 5/50)
g <- as.directed(g)
g <- subgraph(g, subcomponent(g, 1))
graph.cohesion(g)
}
\keyword{graphs}
