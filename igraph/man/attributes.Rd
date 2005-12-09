\name{attributes}
\alias{attributes}
\alias{add.graph.attribute}
\alias{remove.graph.attribute}
\alias{get.graph.attribute}
\alias{set.graph.attribute}
\alias{g.a}
\alias{g.a<-}
\alias{add.vertex.attribute}
\alias{remove.vertex.attribute}
\alias{get.vertex.attribute}
\alias{set.vertex.attribute}
\alias{v.a}
\alias{v.a<-}
\alias{add.edge.attribute}
\alias{remove.edge.attribute}
\alias{get.edge.attribute}
\alias{set.edge.attribute}
\alias{e.a}
\alias{e.a<-}
\title{Graph, vertex and edge attributes}
\description{Attributes associated values belonging to a whole graph,
  its vertices or edges. These can represent some property, like data
  about how the graph was constructed, the color of the vertices when
  the graph is plotted, or simply the weights of the edges in a weighted
  graph.}
\usage{
graph <- add.graph.attribute(graph, attrname)
graph <- remove.graph.attribute(graph, attrname)
get.graph.attribute(graph, attrname=NULL)
set.graph.attribute(graph, attrname, value)
g.a(graph, attrname=NULL)
g.a(graph, attrname) <- value
graph <- add.vertex.attribute(graph, attrname)
graph <- remove.vertex.attribute(graph, attrname)
get.vertex.attribute(graph, attrname=NULL, v=igraph.vs.all(graph))
set.vertex.attribute(graph, attrname, v=igraph.vs.all(graph), value)
v.a(graph, attrname=NULL, v=igraph.vs.all(graph))
v.a(graph, attrname, v=igraph.vs.all(graph)) <- value
graph <- add.edge.attribute(graph, attrname)
graph <- remove.edge.attribute(graph, attrname)
get.edge.attribute(graph, attrname=NULL, e=1:ecount(graph)-1)
set.edge.attribute(graph, attrname, e=1:ecount(graph)-1, value)
e.a(graph, attrname=NULL, e=1:ecount(graph)-1)
e.a(graph, attrname, e=1:ecount(graph)-1) <- value
}
\arguments{
  \item{graph}{The graph object to work on. Note that the original graph
  is never modified, a new graph object is returned instead; if don't
  assign it to a variable your modifications will be lost! See examples
  below.}
  \item{attrname}{Character constant, the name of the attribute. For
    query functions this can be \code{NULL}, in this case the list of
    the attributes will be returned.}
  \item{v}{Numeric vector, the ids of the vertices. It is not recycled,
    even if \code{value} is longer.}
  \item{value}{Numeric vector, the new value(s) of the attributes, it
    will be recycled if needed.}
  \item{e}{Numeric vector, the ids of the edges. It is not recycled,
    even if \code{value} is longer.}
}
\details{Each attribute has a name, this is simply a character
  string. The names of the graph, vertex and edge attributes use
  different name spaces, so the \sQuote{color} vertex attribute has
  nothing to do with the \sQuote{color} edge attribute.

  In order to use an attribute it has to be \sQuote{added} first by
  using the \code{add.*.attribute} functions. The initial values of the
  attributes are undefined. After an attribute is added, its value(s)
  can be set by the \code{set.*.attribute} functions and queried by the
  \code{get.*.attribute} functions.

  You can also use the shorter \code{g.a}, \code{v.a}, \code{e.a} forms,
  see the examples below, they are shorthands for the
  \code{get.*.attribute} and \code{set.*.attribute} functions, the
  latter is in effect if they are used on the left hand side of an
  assignment, see examples below.

  An unneeded attribute can be removed by the \code{remove.*.*}
  functions.

  You can list the attributes of a graph by giving \code{attrname=NULL}
  to the query functions (this is the default), or by
  \code{\link{print.igraph}}.
}
\note{Only numeric constant attributes are implemented right now, string
  attributes will be added shortly.}
\value{
  The \code{add.*.attribute} and \code{remove.*.attribute} and
  \code{set.*.attribute} functions all return a new graph object with the
  requested modifications.

  The \code{get.*.attribute} functions return numeric vectors, the
  values of the attributes if the \code{attrname} argument is not
  \code{NULL} and a character vector, the name of the attributes if it
  is \code{NULL}.
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{print.igraph}} can also print attributes}
\examples{
g <- graph.star(10)
g <- add.graph.attribute(g, "id")
g.a(g, "id")
remove.graph.attribute(g, "id") ## g is unchanged
g.a(g)
g.a(g, "id")

g <- add.vertex.attribute(g, "size")
v.a(g, "size") <- c(1,2)
v.a(g, "size")

g <- add.edge.attribute(g, "weight")
e.a(g, "weight") <- sample(1:100, 20, replace=TRUE)
e.a(g, "weight")

print.igraph(g, g=TRUE, v=TRUE, e=TRUE)

g2 <- graph.star(10)
g2 <- add.graph.attribute(g2, "date", "character")
g.a(g2, "date")
g.a(g2, "date") <- "2005/10/25"
g.a(g2, "date")
remove.graph.attribute(g2, "date") ## g2 is unchanged
g.a(g2)
g.a(g2, "date")

g2 <- add.vertex.attribute(g2, "id", "character")
abc <- strsplit("abcdefghujklmnpqrstuvwxyz", "")[[1]]
v.a(g2, "id") <- abc
v.a(g2, "id")

g2 <- add.edge.attribute(g2, "random", "character")
e.a(g2, "random") <- sample(1:100, 20, replace=TRUE)
e.a(g2, "random")

print.igraph(g2, g=TRUE, v=TRUE, e=TRUE)
}
\keyword{graphs}
