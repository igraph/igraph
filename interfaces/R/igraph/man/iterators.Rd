\name{iterators}
\alias{igraph.es.adj}
\alias{igraph.es.all}
\alias{igraph.es.end}
\alias{igraph.es.from}
\alias{igraph.es.fromorder}
\alias{igraph.es.get}
\alias{igraph.es.next}
\alias{igraph.es.reset}
\alias{$.igraphesseq}
\alias{igraph.es.to}
\alias{$.igraphesvector}
\alias{igraph.es.vector}
\alias{igraph.vs.adj}
\alias{igraph.vs.all}
\alias{igraph.vs.end}
\alias{igraph.vs.get}
\alias{igraph.vs.next}
\alias{igraph.vs.reset}
\alias{$.igraphvsseq}
\alias{$.igraphvsvector}
\alias{igraph.vs.vector}
\alias{as.vector.igraphvsseq}
\alias{as.vector.igraphvsvector}
\alias{as.vector.igraphesseq}
\alias{as.vector.igraphesvector}
\title{Vertex and edge sequences}
\description{Vertex and edge sequences are the central concept of
  igraph.}
\usage{
igraph.vs.all(graph)
igraph.vs.vector(graph, iterator)
igraph.vs.adj(graph, vid, mode="all")

igraph.vs.end(graph, iterator)
igraph.vs.get(graph, iterator)
igraph.vs.next(graph, iterator)
igraph.vs.reset(graph, iterator)

vs$end
vs$step
vs$v

igraph.es.all(graph)
igraph.es.fromorder(graph)
igraph.es.vector(graph, iterator)
igraph.es.adj(graph, vid, mode="all")

igraph.es.end(graph, iterator)
igraph.es.get(graph, iterator)
igraph.es.from(graph, iteraror)
igraph.es.to(graph, iterator)
igraph.es.next(graph, iterator)
igraph.es.reset(graph, iterator)

es$end
es$e
es$from(graph)
es$to(graph)

as.vector(iterator)
}
\arguments{
  \item{graph}{A graph object.}
  \item{iterator}{A vertex or edge sequence.}
  \item{vid}{Id of the vertex of which the adjacent vertices or edges
    will be selected.}
  \item{mode}{Defines the type of adjacent vertices/edges to
    select. \code{out} is for outgoing edges, \code{in} for incoming
    ones, while \code{all} is both both of them.}  
}
\details{
  \code{igraph.vs.all} creates a vertex sequence with all vertices of a
  graph.

  \code{igraph.vs.vector} creates a vertex sequence from a vector of
  vertex ids.

  \code{igraph.vs.adj} creates a vertex sequence containing the adjacent
  vertices of a vertex.

  \code{igraph.vs.end}, \code{igraph.vs.step}, \code{igraph.vs.get}
  provide an iterator abstraction of a vertex set.
  \code{igraph.vs.end} checks whether there is another vertex to step to
  with the iterator. \code{igraph.vs.get} returns the vertex id at the
  current position of the iterator. \code{igraph.vs.step} moves to the
  next position.

  \code{igraph.vs.reset} resets the iterator, this is equivalent to
  creating a new iterator of the same type, with the same parameters.

  \code{vs$end}, \code{vs$step}, \code{vs$v} provide a simpler
  iterator-like interface. See the examples section below.

  \code{as.vector} can be used on a vertex and edge sequence to coerce
  it to a regular numeric vector.

  \code{igraph.es.all} creates an edge sequence containing all edges in
  the graph, their order will be arbitrary.

  \code{igraph.es.fromorder}
  creates an edge sequence with all edges. It is somewhat less efficient
  than \code{igraph.es.all}, but the order of the edges is defined by
  the starting vertex of the edge.

  \code{igraph.es.vector} creates an edge set from a vector containing
  edge ids.

  \code{igraph.es.adj} creates an edge set containing the adjacent edges
  of a vertex.

  \code{igraph.es.get} returns the id of the edge at the current
  iterator position.

  \code{igraph.es.from} returns the starting vertex
  of the current edge.

  \code{igraph.es.to} returns the end vertex of the current edge.

  \code{igraph.es.end} checks whether there are more edges in the edge
  set to step to.

  \code{igraph.es.next} steps to the next edge in the edge sequence.

  \code{igraph.es.reset} resets the iterator, it is equivalent to but slightly
  more efficient than creating a new edge sequence of the same type with
  the same parameters.

  \code{es$end}, \code{es$e}, \code{es$from(graph)} and
  \code{es$to(graph)} provide a simpler iterator interface to edge
  sequences, see the examples below.

  Note that \code{igraph.vs.next}, \code{igraph.vs.reset},
  \code{igraph.es.next} and \code{igraph.es.reset} don't change the
  iterator supplied as a parameter, but return the new, modified
  iterator. 
}
\value{
  \code{igraph.vs.all}, \code{igraph.vs.vector}, \code{igraph.vs.adj}
  return a vertex sequence.

  \code{igraph.vs.end} and \code{igraph.es.end} return a logical
  constant.

  \code{igraph.vs.get}, \code{igraph.es.get}, \code{igraph.es.get},
  \code{igraph.es.from} and \code{igraph.es.to} return a numeric
  constant.

  \code{igraph.vs.next} and \code{igraph.vs.reset} return the new vertex
  sequence.

  \code{igraph.es.next} and \code{igraph.es.reset} return the new vertex
  sequence.

  \code{as.vector} returns a numeric vector.
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
# mean degree of vertices in the largest cluster in a random graph
g <- erdos.renyi.game(100, 2/100)
c <- clusters(g)
vsl <- igraph.vs.vector(g, which(which.max(c$csize)-1==c$membership)-1)
mean(degree(g, vsl))

# set the color of these vertices to red, others greens
g <- add.vertex.attribute(g, "color", type="string")
v.a(g, "color") <- "green"
v.a(g, "color", vsl) <- "red"
\dontrun{plot(g, vertex.size=3, labels=NA, vertex.color="a:color")}

# the longest geodesic within the largest cluster
vsl <- igraph.vs.reset(g, vsl) # just to be sure
long <- numeric()
while (!vsl$end) {
  paths <- get.shortest.paths(g, vsl$v)
  paths <- paths[ as.vector(vsl)+1  ]
  fl <- paths[[ which.max(sapply(paths, length)) ]]
  if (length(fl) > length(long)) {
    long <- fl
  }
  vsl <- vsl$step
}
}
\keyword{graphs}
