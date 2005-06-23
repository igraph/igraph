\name{attributes}
\alias{add.graph.attribute}
\alias{delete.graph.attribute}
\alias{get.graph.attribute}
\alias{set.graph.attribute}
\alias{add.vertex.attribute}
\alias{delete.vertex.attribute}
\alias{get.vertex.attribute}
\alias{set.vertex.attribute}
\alias{add.edge.attribute}
\alias{delete.edge.attribute}
\alias{get.edge.attribute}
\alias{set.edge.attribute}
\alias{g.a}
\alias{v.a}
\alias{e.a}
\title{Graph, vertex and edge attributes}
\description{Graph attributes are arbitrary valued assigned to graph
  objects, while vertex attributes are assigned to vertices and edge
  attributes to edges.}
\usage{
add.graph.attribute(graph, attrname, default = NA)
delete.graph.attribute(graph, attrname)
get.graph.attribute(graph, attrname)
set.graph.attribute(graph, attrname, value)
add.vertex.attribute(graph, attrname, default = NA)
delete.vertex.attribute(graph, attrname)
get.vertex.attribute(graph, attrname, v = NULL)
set.vertex.attribute(graph, attrname, v = NULL, value)
add.edge.attribute(graph, attrname, default = NA)
delete.edge.attribute(graph, attrname)
get.edge.attribute(graph, attrname, from = NULL, to = NULL)
set.edge.attribute(graph, attrname, from = NULL, to = NULL, value)
g.a(graph, attrname)
g.a(graph, attrname) <- value
v.a(graph, attrname, v = NULL)
v.a(graph, attrname, v = NULL) <- value
e.a(graph, attrname, from = NULL, to = NULL)
e.a(graph, attrname, from = NULL, to = NULL) <- value
}
\arguments{
  \item{graph}{The graph to work on.}
  \item{attrname}{Character constant, the name of the attribute.}
  \item{default}{The default value of the attribute, arbitrary object.}
  \item{value}{The new value of the attribute.}
  \item{from}{The tail of the edge.}
  \item{to}{The head of the edge.}
}
\details{This manual page is not yet (well) written, as the attribute
  handling functions, especially those for edge attributes, are
  expected to be rewised shortly.

  If order to assign (any type of) attributes, first you have to call
  \code{add.*.attribute} and then \code{set.*.attribute} (or one of its
  shorter forms).

  Adding a vertex or edge attribute to a graph is a cheap
  operation. Setting attributes however might be costly.
}
\value{The value(s) of the attributes for query functions, an updated
  graph object for the others.}
\note{Edge attributes are not very well designed to use with graphs
  containing multiple edges.}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
g <- graph.ring(10, name="ring")
g <- add.vertex.attribute(g, "color", default="red")
g <- set.vertex.attribute(g, "color", 1:5, "green")
v.a(g, "color")
}
\keyword{graphs}
