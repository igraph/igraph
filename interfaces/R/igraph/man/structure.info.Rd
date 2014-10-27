\name{structure.info}
\alias{gorder}
\alias{vcount}
\alias{gsize}
\alias{ecount}
\alias{neighbors}
\alias{incident}
\alias{is.directed}
\alias{is_directed}
\alias{is_adjacent_to}
\alias{are.connected}
\alias{get.edge}
\alias{get.edges}
\alias{ends}
\alias{head_of}
\alias{tail_of}
\title{Gaining information about graph structure}
\description{Functions for exploring the basic structure of a network:
  number of vertices and edges, the neighbors of a node, test whether
  two vertices are connected by an edge.
}
\usage{
gorder(graph)
vcount(graph)
gsize(graph)
ecount(graph)
neighbors(graph, v, mode = 1)
incident(graph, v, mode=c("all", "out", "in", "total"))
is_directed(graph)
is_adjacent_to(graph, v1, v2)
get.edge(graph, id)
ends(graph, es, names = TRUE)
head_of(graph, es, names = TRUE)
tail_of(graph, es, names = TRUE)
}
\arguments{
  \item{graph}{The graph.}
  \item{v}{The vertex of which the adjacent vertices or incident edges
    are queried.}
  \item{mode}{Character string, specifying the type of adjacent vertices
    or incident edges to list in a directed graph. If \dQuote{out}, then
    only outgoing edges (or their corresponding vertices) are
    considered; \dQuote{in} considers incoming edges; \sQuote{all}
    ignores edge directions. This argument is ignored for undirected
    graphs.}
  \item{v1}{The id of the first vertex. For directed graphs only edges
    pointing from \code{v1} to \code{v2} are searched.}
  \item{v2}{The id of the second vertex. For directed graphs only edges
    pointing from \code{v1} to \code{v2} are searched.}
  \item{id}{A numeric edge id.}
  \item{es}{An edge sequence.}
  \item{names}{Whether to return symbolic names for named graphs.}
}
\details{
  These functions provide the basic structural information of a graph.

  \code{gorder} gives the order of the graph, the number of vertices.
  \code{vcount} is an alias of this function. \code{gorder} is available
  from igraph version 0.8.0.

  \code{gsize} gives the size of the graph, its number of edges.
  \code{ecount} is an alias of \code{gsize}, and is available from
  igraph version 0.8.0.

  \code{neighbors} gives the neighbors of a vertex. The vertices
  connected by multiple edges are listed as many times as the number of
  connecting edges.

  \code{incident} gives the incident edges of a vertex.
  
  \code{is_directed} gives whether the graph is directed or not. It just
  gives its \code{directed} attribute.
  
  \code{is_adjacent_to} decides whether there is an edge from \code{v1}
  to \code{v2}. 

  \code{get.edge} returns the end points of the edge with the supplied
  edge id. For directed graph the source vertex comes first, for
  undirected graphs, the order is arbitrary. Note that \code{get.edge}
  is deprecated from igraph version 0.8.0 and will be eventually
  removed from igraph. Please use \code{ends} instead.

  \code{ends} returns a matrix with the endpoints of the edges in
  the edge sequence argument.

  \code{head_of} returns the head vertices of directed edges,
  \code{tail_of} returns their tail vertices. They return vertex
  sequences. (For undirected graphs, they return an arbitrary incident
  vertex for each edge, but together they return both incident vertices
  for each edge.)

}
\value{
  \code{gorder} and \code{gsize} return integer
  constants. \code{neighbors} returns a vertex sequence.
  \code{is_directed} and \code{is_adjacent_to} return boolean
  constants. \code{get.edge} returns a numeric vector of length two.
  \code{ends} returns a two-column numeric matrix.
}
%\references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{graph}}}
\examples{
g <- make_ring(10)
gorder(g)
vcount(g)
gsize(g)
ecount(g)
neighbors(g, 5)
incident(g, 5)
is_adjacent_to(g, 1, 2)
is_adjacent_to(g, 2, 4)
ends(g, 1:6)

g2 <- make_directed_graph(c("A", "B"))
head_of(g2, 1)
tail_of(g2, 1)
}
\keyword{graphs}
