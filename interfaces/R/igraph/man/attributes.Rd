\name{attributes}
\alias{attributes}
\alias{set.graph.attribute}
\alias{set_graph_attr}
\alias{get.graph.attribute}
\alias{graph_attr}
\alias{graph_attr<-}
\alias{remove.graph.attribute}
\alias{delete_graph_attr}
\alias{list.graph.attributes}
\alias{graph_attr_names}
\alias{graph.attributes}
\alias{graph.attributes<-}
\alias{set.vertex.attribute}
\alias{set_vertex_attr}
\alias{get.vertex.attribute}
\alias{vertex_attr}
\alias{vertex_attr<-}
\alias{remove.vertex.attribute}
\alias{delete_vertex_attr}
\alias{list.vertex.attributes}
\alias{vertex_attr_names}
\alias{vertex.attributes}
\alias{vertex.attributes<-}
\alias{set.edge.attribute}
\alias{set_edge_attr}
\alias{get.edge.attribute}
\alias{edge_attr}
\alias{edge_attr<-}
\alias{remove.edge.attribute}
\alias{delete_edge_attr}
\alias{list.edge.attributes}
\alias{edge_attr_names}
\alias{edge.attributes}
\alias{edge.attributes<-}
\concept{Vertex/edge/graph attributes}
\title{Graph, vertex and edge attributes}
\description{Attributes are associated values belonging to a graph,
  vertices or edges. These can represent some property, like data
  about how the graph was constructed, the color of the vertices when
  the graph is plotted, or simply the weights of the edges in a weighted
  graph.}
\usage{
graph_attr(graph, name)
set_graph_attr(graph, name, value)
graph_attr_names(graph)
graph.attributes(graph)
delete_graph_attr(graph, name)
vertex_attr(graph, name, index=V(graph))
set_vertex_attr(graph, name, index=V(graph), value)
delete_vertex_attr(graph, name)
vertex_attr_names(graph)
vertex.attributes(graph, index = V(graph))
edge_attr(graph, name, index=E(graph))
set_edge_attr(graph, name, index=E(graph), value)
delete_edge_attr(graph, name)
edge_attr_names(graph)
edge.attributes(graph, index = E(graph))
}
\arguments{
  \item{graph}{The graph object to work on. Note that the original graph
  is never modified, a new graph object is returned instead; if you don't
  assign it to a variable your modifications will be lost! See examples
  below.}
  \item{name}{Character constant, the name of the attribute.}
  \item{index}{Numeric vector, the ids of the vertices or edges.
    It is not recycled, even if \code{value} is longer.}
  \item{value}{Numeric vector, the new value(s) of the attributes, it
    will be recycled if needed.}
}
\details{
  There are three types of attributes in igraph: graph, vertex and edge
  attributes. Graph attributes are associated with graph, vertex
  attributes with vertices and edge attributes with edges.

  Examples for graph attributes are the date when the graph data was
  collected or other types of memos like the type of the data, or
  whether the graph is a simple graph, ie. one without loops and
  multiple edges.

  Examples of vertex attributes are vertex properties, like the
  vertex coordinates for the visualization of the graph, or other
  visualization parameters, or meta-data associated with the vertices,
  like the gender and the age of the individuals in a friendship
  network, the type of the neurons in a graph representing neural
  circuitry or even some pre-computed structual properties, like the
  betweenness centrality of the vertices.

  Examples of edge attributes are data associated with edges: most
  commonly edge weights, or visualization parameters.

  In recent igraph versions, arbitrary R objects can be assigned as
  graph, vertex or edge attributes.
  
  Some igraph functions use the values or graph, vertex and edge
  attributes if they are present but this is not done in the current
  version very extensively. Expect more in the (near) future.

  Graph attributes can be created with the \code{set_graph_attr}
  function, and removed with \code{delete_graph_attr}. Graph
  attributes are queried with \code{graph_attr} and the
  assigned graph attributes are listed with
  \code{graph_attr_names}.

  If no attribute name is specified, \code{graph_attr} returns all
  graph attributes. This function also has a counterpart that sets all
  graph attributes at once, see an example below.

  There is a simpler notation for using graph attributes: the
  \sQuote{\code{$}} operator. It can be used both to query and set graph
  attributes, see an example below.
  
  The functions for vertex attributes are
  \code{set_vertex_attr}, \code{vertex_attr},
  \code{delete_vertex_attr} and \code{vertex_attr_names} and
  for edge attributes they are \code{set_edge_attr},
  \code{edge_attr}, \code{delete_edge_attr} and
  \code{edge_attr_names}.

  Similarly to graph attributes, \code{vert_attr} returns all
  vertex attributes if no attribute name is given, in a named list, and
  \code{edge_attr} returns all edge attributes, in a named list.
  
  There is however a (syntactically) much simpler way to handle vertex
  and edge attribute by using vertex and edge selectors, it works like this:
  \code{V(g)} selects all vertices in a graph, and \code{V(g)$name}
  queries the \code{name} attribute for all vertices. Similarly is
  \code{vs} is a vertex set \code{vs$name} gives the values of the
  \code{name} attribute for the vertices in the vertex set.

  This form can also be used to set the values of the attributes, like
  the regular R convention:
  \preformatted{V(g)$color <- "red"}
  It works for vertex subsets as well:
  \preformatted{V(g)[1:5]$color <- "green"}

  The notation for edges is similar: \code{E(g)} means all edges
  \code{E(g)$weight} is the \code{weight} attribute for all edges, etc.

  See also the manual page for \code{iterators} about how to create
  various vertex and edge sets.
}
\value{
  \code{graph_attr}, \code{vertex_attr} and
  \code{edge_attr} return an R object, or a list of R objects
  if attributes of more vertices/edges are requested.

  \code{set_graph_attr}, \code{set_vertex_attr},
  \code{set_edge_attr}, and also \code{delete_graph_attr},
  \code{delete_vertex_attr} and \code{delete_edge_attr} return
  a new graph object with the updates/removes performed.

  \code{graph_attr_names}, \code{vertex_attr_names} and
  \code{edge_attr_names} return a character vector, the names of
  the attributes present.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{print.igraph}} can print attributes. See
  \code{\link{attribute.combination}} for details on how igraph combines
  attributes if several vertices or edges are mapped into one.}
\examples{
g <- make_ring(10)
g <- set_graph_attr(g, "name", "RING")
# It is the same as
g$name <- "RING"
g$name

g <- set_vertex_attr(g, "color", value=c("red", "green"))
vertex_attr(g, "color")

g <- set_edge_attr(g, "weight", value=runif(ecount(g)))
edge_attr(g, "weight")

# The following notation is more convenient
g <- make_star(10)

V(g)$color <- c("red", "green")
V(g)$color

E(g)$weight <- runif(ecount(g))
E(g)$weight

str(g, g=TRUE, v=TRUE, e=TRUE)

# Setting all attributes at once
g2 <- make_empty_graph(10)
g2

graph_attr(g2) <- graph_attr(g)
vertex_attr(g2) <- vertex_attr(g)
edge_attr(g2) <- list()
g2

edge_attr(g2) <- list(weight=numeric())
g2
}
\keyword{graphs}
