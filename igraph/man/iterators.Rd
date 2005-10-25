\name{iterators}
\alias{iterators}
\alias{ii.create}
\alias{ii.end}
\alias{ii.get.edge}
\alias{ii.get.from}
\alias{ii.get.to}
\alias{ii.get.vertex}
\alias{ii.get.vertex.nei}
\alias{ii.next}
\alias{ii.prev}
\alias{ii.reset}
\title{Vertex and edge iterators}
\description{Iterators provide a handy method to walk through the edges
  or vertices of a graph.}
\usage{
iterator <- ii.create(graph, type, vid=NULL, mode=NULL)
iterator <- ii.next(graph, iterator)
iterator <- ii.prev(graph, iterator)
iterator <- ii.reset(graph, iterator)
ii.end(graph, iterator)
ii.get.vertex(graph, iterator)
ii.get.edge(graph, iterator)
ii.get.from(graph, iterator)
ii.get.to(graph, iterator)
ii.get.vertex.nei(graph, iterator)
}
\arguments{
  \item{graph}{The graph object to walk through. Note that iterator
    object does not store its graph, or a reference to it, so this
    argument has to be given for every function documented here. }
  \item{type}{The type of the iterator, see details below.}
  \item{vid}{The vertex id for the \sQuote{eneis} and \sQuote{vneis}
    iterator types, see details below.}
  \item{mode}{The mode for the \sQuote{eneis} and \sQuote{vneis}
    iterator types, see details below.}
  \item{iterator}{The iterator object.}
}
\details{
  There are two types of iterators: vertex iterators walk through
  vertices of a graph, edge iterators walk through edges. Vertex
  iterators does not neccesserily walk through all edges, and similarly
  edge iterators might use only a subset of all edges in the graph.  
  
  To create an iterator, one should call \code{ii.create} first, giving
  the graph an the type of the iterator. See the specific iterator types
  below. If the structure of the graph
  is modified after creating the iterator, the iterator becomes invalid
  and the result of calling any function on it is undefined!

  After creation, the iterator points to the first vertex/edge in the
  graph, and \code{ii.next} moves it to the next
  vertex/edge. \code{ii.end} checks whether there are more
  vertices/edges to visit. For vertex iterators \code{ii.get.vertex}
  gives the id of the current vertex. For edge iterators
  \code{ii.get.edge} gives the current edge and \code{ii.get.from} and
  \code{ii.get.to} give the two end points of the current edge. For
  directed graphs \code{ii.get.from} always gives the starting point and
  \code{ii.get.to} always gives the end point of the edge.

  Vertex iterators:
  \itemize{
    \item{\sQuote{vid},}{simple vertex iterator, steps through the
      vertices according to increasing vertex ids.}
    \item{\sQuote{vneis},}{steps through the neighboring vertices of a
      given vertex (the \code{vid} argument). In directed graphs the
      \code{mode} argument specifies the neighbors to step through, see
      the same argument of the \code{\link{neighbors}} function for
      details.}
  }

  Edge iterators:
  \itemize{
    \item{\sQuote{eid},}{simple edge iterator, steps through the edges
      according to increasing edge ids.}
    \item{\sQuote{efromorder},}{steps through edges according to the
      increasing vertex id of the starting point of the edges. This only
      makes sense in directed graphs, for undirected graphs, the order
      will be arbitrary.}
    \item{\sQuote{eneis},}{steps through the edges attached to a
      vertex. The vertex is given by the \code{vid} argument, while the
      \code{mode} argument specifies the type of edges (incoming and/or
      outgoing) to consider. See the \code{mode} argument of the
      \code{\link{neighbors}} function for details.}
  }

  Some other rarely used iterator functions: \code{ii.reset} can reset
  an iterator to its first vertex/edge. \code{ii.prev} goes to the
  previous vertex/edge. \code{ii.get.vertex.nei} gives the id of the
  \dQuote{other} vertex at the end of the edge for the \sQuote{eneis}
  iterator.
}
\note{Not all functions are implemented for all iterators. Some make
  sense only for vertex or edge iterators, other (like \code{ii.prev})
  are rarely used, and rarely implemented.
}
\value{
  \code{ii.create} returns an iterator object.

  \code{ii.next}, \code{ii.prev}, \code{ii.reset} return a new, modified
  (stepped, reseted, etc.) iterator object.

  \code{ii.end} returns a logical constant.

  \code{ii.get.vertex}, \code{ii.get.from}, \code{ii.get.to} and
  \code{ii.get.vertex.nei} return a vertex id, a numeric constant.

  \code{ii.get.edge} returns an edge id, a numeric constant.
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
# iterate through the edges of a graph
g <- graph.star(10)
it <- ii.create(g, "eid")
while (!ii.end(g, it)) {
  print(paste( ii.get.from(g, it), ii.get.to(g, it) ))
  it <- ii.next(g, it)
}
}
\keyword{graphs}
