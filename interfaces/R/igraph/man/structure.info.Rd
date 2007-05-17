\name{structure.info}
\alias{vcount}
\alias{ecount}
\alias{neighbors}
\alias{is.directed}
\alias{are.connected}
\alias{get.edge}
\alias{get.edges}
\title{Gaining information about graph structure}
\description{Functions for exploring the basic structure of a network:
  number of vertices and edges, the neighbors of a node, test whether
  two vertices are connected by an edge.
}
\usage{
vcount(graph)
ecount(graph)
neighbors(graph, v, mode = 1)
is.directed(graph)
are.connected(graph, v1, v2)
get.edge(graph, id)
get.edges(graph, es)
}
\arguments{
  \item{graph}{The graph.}
  \item{v}{The vertex of which the neighbors are queried.}
  \item{mode}{Character string, specifying the type of neighboring
    vertices to list in a directed graph. If \dQuote{out} the vertices
    \emph{to} which an edge exist are listed, if \dQuote{in} the
    vertices \emph{from} which an edge is directed are listed,
    \dQuote{all} lists all the vertices. This argument is ignored for
    undirected graphs.} 
  \item{v1}{The id of the first vertex. For directed graphs only edges
    pointing from \code{v1} to \code{v2} are searched.}
  \item{v2}{The id of the second vertex. For directed graphs only edges
    pointing from \code{v1} to \code{v2} are searched.}
  \item{id}{A numeric edge id.}
  \item{es}{An edge sequence.}
}
\details{
  These functions provide the basic structural information of a graph.

  \code{vcount} gives the number of vertices in the graph.

  \code{ecount} gives the number of edges in the graph.

  \code{neighbors} gives the neighbors of a vertex. The vertices
  connected by multiple edges are listed as many times as the number of
  connecting edges.
  
  \code{is.directed} gives whether the graph is directed or not. It just
  gives its \code{directed} attribute.
  
  \code{are.connected} decides whether there is an edge from \code{v1}
  to \code{v2}. 

  \code{get.edge} returns the end points of the edge with the supplied
  edge id. For directed graph the source vertex comes first, for
  undirected graphs, the order is arbitrary.

  \code{get.edges} returns a matrix with the endpoints of the edges in
  the edge sequence argument. 
}
\value{
  \code{vcount} and \code{ecount} return integer
  constants. \code{neighbors} returns an integer
  vector. \code{is.directed} and \code{are.connected} return boolean
  constants. \code{get.edge} returns a numeric vector of length two.
  \code{get.edges} returns a two-column matrix.
}
%\references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{graph}}}
\examples{
g <- graph.ring(10)
vcount(g)
ecount(g)
neighbors(g, 5)
are.connected(g, 1, 2)
are.connected(g, 2, 4)
get.edges(g, 0:5)
}
\keyword{graphs}
