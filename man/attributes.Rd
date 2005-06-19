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
\description{}
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
  \item{graph}{}
  \item{attrname}{}
  \item{default}{}
  \item{value}{}
  \item{from}{}
  \item{to}{}
}
\details{}
\value{}
\references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{}
\examples{}
\keyword{graphs}
