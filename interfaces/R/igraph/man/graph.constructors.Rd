\name{graph.constructors}
\alias{graph.constructors}
\alias{graph.empty}
\alias{graph}
\alias{graph.star}
\alias{graph.lattice}
\alias{graph.ring}
\alias{graph.tree}
\alias{graph.full}
\alias{graph.full.citation}
\alias{graph.atlas}
\alias{graph.edgelist}
\alias{graph.extended.chordal.ring}
\concept{Tree}
\concept{Lattice}
\concept{Star graph}
\concept{Graph Atlas}
\concept{Empty graph}
\concept{Full graph}
\title{Various methods for creating graphs}
\description{These method can create various (mostly regular) graphs:
  empty graphs, graphs with the given edges, graphs from adjacency
  matrices, star graphs, lattices, rings, trees.}
\synopsis{
graph.lattice(dimvector = NULL, length = NULL, dim = NULL, nei = 1,
              directed = FALSE, mutual = FALSE, circular = FALSE, \dots)
}
\usage{
graph.empty(n=0, directed=TRUE)
graph(edges, n=max(edges)+1, directed=TRUE)
graph.star(n, mode = "in", center = 0)
graph.lattice(dimvector, nei = 1, directed = FALSE, mutual = FALSE, 
              circular = FALSE)
graph.lattice(length, dim, nei = 1, directed = FALSE, mutual = FALSE,
              circular = FALSE)
graph.ring(n, directed = FALSE, mutual = FALSE, circular=TRUE)
graph.tree(n, children = 2, mode="out")
graph.full(n, directed = FALSE, loops = FALSE)
graph.full.citation(n, directed = TRUE)
graph.atlas(n)
graph.edgelist(el, directed=TRUE)
graph.extended.chordal.ring(n, w)
}
\arguments{
  \item{edges}{Numeric vector defining the edges, the first edge points
    from the first element to the second, the second edge from the third
  to the fourth, etc.}
  \item{directed}{Logical, if TRUE a directed graph will be
    created. Note that for while most constructors the default is TRUE,
    for \code{graph.lattice} and \code{graph.ring} it is FALSE. For
    \code{graph.star} the \code{mode} argument should be used for
    creating an undirected graph.}
  \item{n}{The number of vertices in the graph for most functions.

    For \code{graph} this parameter is ignored if there is a bigger
    vertex id in \code{edges}. This means that for this function it is
    safe to supply zero here if the vertex with the largest id is not an
    isolate.

    For \code{graph.atlas} this is the number (id) of the graph to
    create. 
  }
  \item{mode}{
    For \code{graph.star} it defines the direction of the
    edges, \code{in}: the edges point \emph{to} the center, \code{out}:
    the edges point \emph{from} the center, \code{undirected}: the edges
    are undirected.

    For \code{igraph.tree} this parameter defines the direction of the
    edges. \code{out} indicates that the edges point from the parent to
    the children, \code{in} indicates that they point from the children
    to their parents, while \code{undirected} creates an undirected
    graph. 
  }
  \item{center}{For \code{graph.star} the center vertex of the graph, by
    default the first vertex.}
  \item{dimvector}{A vector giving the size of the lattice in each
    dimension, for \code{graph.lattice}.}
  \item{nei}{The distance within which (inclusive) the neighbors on the
    lattice will be connected. This parameter is not used right now.}
  \item{mutual}{Logical, if TRUE directed lattices will be mutually
    connected.} 
  \item{circular}{Logical, if TRUE the lattice or ring will be circular.}
  \item{length}{Integer constant, for regular lattices, the size of the
    lattice in each dimension.}
  \item{dim}{Integer constant, the dimension of the lattice.}
  \item{children}{Integer constant, the number of children of a vertex
    (except for leafs) for \code{graph.tree}.}
  \item{loops}{If TRUE also loops edges (self edges) are added.}
  \item{graph}{An object.}
  \item{el}{An edge list, a two column matrix, character or numeric. See
    details below.}
  \item{w}{A matrix which specifies the extended chordal ring. See
    details below.}
}
\details{All these functions create graphs in a deterministic way.

  \code{graph.empty} is the simplest one, this creates an empty graph.

  \code{graph} creates a graph with the given edges.

  \code{graph.star} creates a star graph, in this every single vertex is
  connected to the center vertex and nobody else.

  \code{graph.lattice} is a flexible function, it can create lattices of
  arbitrary dimensions, periodic or unperiodic ones.

  \code{graph.ring} is actually a special case of \code{graph.lattice},
  it creates a one dimensional circular lattice.

  \code{graph.tree} creates regular trees.

  \code{graph.full} simply creates full graphs.

  \code{graph.full.citation} creates a full citation graph. This is a
  directed graph, where every i->j edge is present if and only if j<i.
  If \code{directed=FALSE} then the graph is just a full graph.
  
  \code{graph.atlas} creates graphs from the book An Atlas of Graphs by
  Roland C. Read and Robin J. Wilson. The atlas contains all undirected
  graphs with up to seven vertices, numbered from 0 up to 1252. The
  graphs are listed: 
  \enumerate{
    \item in increasing order of number of nodes;
    \item for a fixed number of nodes, in increasing order of the number
      of edges;
    \item for fixed numbers of nodes and edges, in increasing order of
      the degree sequence, for example 111223 < 112222;
    \item for fixed degree sequence, in increasing number of
      automorphisms.
  }

  \code{graph.edgelist} creates a graph from an edge list. Its argument
  is a two-column matrix, each row defines one edge. If it is 
  a numeric matrix then its elements are interpreted as vertex ids. If
  it is a character matrix then it is interpreted as symbolic vertex
  names and a vertex id will be assigned to each name, and also a
  \code{name} vertex attribute will be added.
  
  \code{graph.extended.chordal.ring} creates an extended chordal ring.
  An extended chordal ring is regular graph, each node has the same
  degree. It can be obtained from a simple ring by adding some extra
  edges specified by a matrix. Let p denote the number of columns in
  the \sQuote{\code{W}} matrix. The extra edges of vertex \code{i}
  are added according to column \code{i mod p} in
  \sQuote{\code{W}}. The number of extra edges is the number
  of rows in \sQuote{\code{W}}: for each row \code{j} an edge
  \code{i->i+w[ij]} is added if \code{i+w[ij]} is less than the number
  of total nodes. See also Kotsis, G: Interconnection Topologies for
  Parallel Processing Systems, PARS Mitteilungen 11, 1-6, 1993.

  
}
\value{Every function documented here returns a \code{graph} object.}
%\references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{graph.adjacency}} to create graphs from adjacency
  matrices, \code{\link{graph.formula}} for a handy way to create small
  graphs, \code{\link{graph.data.frame}} for an easy way to create
  graphs with many edge/vertex attributes.
}
\examples{
g1 <- graph.empty()
g2 <- graph( c(1,2,2,3,3,4,5,6), directed=FALSE )
g5 <- graph.star(10, mode="out")
g6 <- graph.lattice(c(5,5,5))
g7 <- graph.lattice(length=5, dim=3)
g8 <- graph.ring(10)
g9 <- graph.tree(10, 2)
g10 <- graph.full(5, loops=TRUE)
g11 <- graph.full.citation(10)
g12 <- graph.atlas(sample(0:1252, 1))
el <- matrix( c("foo", "bar", "bar", "foobar"), ncol=2, byrow=TRUE)
g13 <- graph.edgelist(el)
g15 <- graph.extended.chordal.ring(15, matrix(c(3,12,4,7,8,11), nrow=2))
}
\keyword{graphs}
