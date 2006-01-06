\name{graph.constructors}
\alias{graph.empty}
\alias{graph}
\alias{graph.adjacency}
\alias{graph.star}
\alias{graph.lattice}
\alias{graph.ring}
\alias{graph.tree}
\alias{graph.full}
\title{Various methods for creating graphs}
\description{These method can create various (mostly regular) graphs:
  empty graphs, graphs with the given edges, graphs from adjacency
  matrices, star graphs, lattices, rings, trees.}
\synopsis{
graph.lattice(dimvector = NULL, length = NULL, dim = NULL, nei = 1,
              directed = FALSE, mutual = FALSE, circular = FALSE, \dots)
}
\usage{
graph.empty(\dots, type = igraph.par("default.type"))
graph(edges, \dots)
graph.adjacency(adjmatrix, directed = TRUE, \dots)
graph.star(n, mode = "in", center = 1, directed = TRUE, \dots)
graph.lattice(dimvector, nei = 1, directed = FALSE, mutual = FALSE, 
              circular = FALSE, \dots)
graph.lattice(length, dim, nei = 1, directed = FALSE, mutual = FALSE,
              circular = FALSE, \dots)
graph.ring(n, circular = TRUE, directed = FALSE, mutual = FALSE, \dots)
graph.tree(n, children = 2, \dots)
graph.full(n, directed = FALSE, loops = FALSE)
}
\arguments{
  \item{type}{The storage type of the graph, this defines how its data
    will be stored in the memory or on the disk. Currently the following
    storage types are defined:
    \enumerate{
      \item{\code{adjacencylist}}{the neighbor of each
	node is stored in a list, this type is immutable, graph
	modifying operations create a new graph}
    }
  }
  \item{edges}{Numeric vector defining the edges, the first edge points
    from the first element to the second, the second edge from the third
  to the fourth, etc.}
  \item{adjmatrix}{A square adjacency matrix.}
  \item{directed}{Logical, if TRUE a directed graph will be
    created. Note that for while most constructors the default is TRUE,
    for \code{graph.lattice} and \code{graph.ring} it is FALSE. For
    \code{graph.star} the \code{mode} argument should be used for
    creating an undirected graph.}
  \item{n}{The number of vertices in the graph.}
  \item{mode}{For \code{graph.star} it defined the direction of the
    edges, \code{in}: the edges point \emph{to} the center, \code{out}:
    the edges point \emph{from} the center, \code{undirected}: the edges
    are undirected.}
  \item{center}{For \code{graph.star} the center vertex of the graph, by
    default the first vertex.}
  \item{dimvector}{A vector giving the size of the lattice in each
    dimension, for \code{graph.lattice}.}
  \item{nei}{The distance within which (inclusive) the neighbors on the
    lattice will be connected.}
  \item{mutual}{Logical, if TRUE directed lattices we be mutually
    connected.} 
  \item{circular}{Logical, if TRUE the lattice will be circular.}
  \item{length}{Integer constant, for regular lattices, the size of the
    lattice in each dimension.}
  \item{dim}{Integer constant, the dimension of the lattice.}
  \item{children}{Integer constant, the number of children of a vertex
    for \code{graph.tree}.}
  \item{loops}{If TRUE also loops edges (self edges) are added.}
  \item{graph}{An object.}
  \item{\dots}{Additional arguments, ignored right now.}
}
\details{All these functions create graphs in a deterministic way.

  \code{graph.empty} is the simplest one, this creates an empty graph.

  \code{graph} creates a graph with the given edges.

  \code{graph.adjacency} creates a graph from an adjacency matrix.

  \code{graph.star} creates a star graph, in this every single vertex is
  connected to the center vertex and nobody else.

  \code{graph.lattice} is a flexible function, it can create lattices of
  arbitrary dimensions, periodic or unperiodic ones.

  \code{graph.ring} is actually a special case of \code{graph.lattice},
  it creates a one dimensional circular lattice.

  \code{graph.tree} creates regular trees.

  \code{graph.full} simply creates full graphs.

}
\value{Every function documented here returns a \code{graph} object.}
%\references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
%\seealso{}
\examples{
g1 <- graph.empty()
g2 <- graph( c(1,2,2,3,3,4,5,6), directed=FALSE )
adjm <- matrix(sample(0:1, 100, replace=TRUE, prob=c(0.9,0.1)), nc=10)
g3 <- graph.adjacency( adjm )
g4 <- graph.star(10, mode="out")
g5 <- graph.lattice(c(5,5,5))
g6 <- graph.lattice(length=5, dim=3)
g7 <- graph.ring(10)
g8 <- graph.tree(10, 2)
g9 <- graph.full(5, loops=TRUE)
}
\keyword{graphs}
