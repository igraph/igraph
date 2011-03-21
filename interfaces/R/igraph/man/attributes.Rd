\name{attributes}
\alias{attributes}
\alias{set.graph.attribute}
\alias{get.graph.attribute}
\alias{remove.graph.attribute}
\alias{list.graph.attributes}
\alias{set.vertex.attribute}
\alias{get.vertex.attribute}
\alias{remove.vertex.attribute}
\alias{list.vertex.attributes}
\alias{set.edge.attribute}
\alias{get.edge.attribute}
\alias{remove.edge.attribute}
\alias{list.edge.attributes}
\concept{Vertex/edge/graph attributes}
\title{Graph, vertex and edge attributes}
\description{Attributes are associated values belonging to a graph,
  vertices or edges. These can represent some property, like data
  about how the graph was constructed, the color of the vertices when
  the graph is plotted, or simply the weights of the edges in a weighted
  graph.}
\usage{
get.graph.attribute(graph, name)
set.graph.attribute(graph, name, value)
remove.graph.attribute(graph, name)
get.vertex.attribute(graph, name, index=V(graph))
set.vertex.attribute(graph, name, index=V(graph), value)
remove.vertex.attribute(graph, name)
get.edge.attribute(graph, name, index=E(graph))
set.edge.attribute(graph, name, index=E(graph), value)
remove.edge.attribute(graph, name)
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

  Graph attributes can be created with the \code{set.graph.attribute}
  function, and removed with \code{remove.graph.attribute}. Graph
  attributes are queried with \code{get.graph.attribute} and the
  assigned graph attributes are listed with
  \code{list.graph.attributes}.

  Note that \code{set.graph.attribute} and
  \code{remove.graph.attribute} both return a new graph, they do
  \emph{not} change their argument! If you don't want to lose the
  results, you have to assign their return value to a
  variable: \preformatted{  G <- set.graph.attribute(G, "name", "My favorite graph")
  G <- remove.graph.attribute(G, "name")}
  There is a simpler notation for using graph attributes: the
  \sQuote{\code{$}} operator. It can be used both to query and set graph
  attributes, see an example below.
  
  The functions for vertex attributes are
  \code{set.vertex.attribute}, \code{get.vertex.attribute},
  \code{remove.vertex.attribute} and \code{list.vertex.attributes} and
  for edge attributes they are \code{set.edge.attribute},
  \code{get.edge.attribute}, \code{remove.edge.attribute} and
  \code{list.edge.attributes}.

  Note that \code{set.vertex.attribute},
  \code{remove.vertex.attribute}, \code{set.edge.attribute} and
  \code{remove.edge.attribute} do \emph{not} change their argument,
  you have to assign their return value to a
  variable: \preformatted{  G <- set.vertex.attribute(G, "label", value=letters[1:10])
  G <- remove.vertex.attribute(G, "label")
  G <- set.edge.attribute(G, "weight", value=runif(ecount(G)))
  G <- remove.edge.attribute(G, "weight")}
  
  There is however a (syntactically) much simpler way to handle vertex
  and edge attribute by using vertex and edge selectors, it works like this:
  \code{V(g)} selects all vertices in a graph, and \code{V(g)$name}
  queries the \code{name} attribute for all vertices. Similarly is
  \code{vs} is a vertex set \code{vs$name} gives the values of the
  \code{name} attribute for the vertices in the vertex set.

  This form can also be used to set the values of the attributes, like
  the regular R
  convention: \preformatted{V(g)$color <- "red"}
  It works for vertex subsets as
  well: \preformatted{V(g)[0:5]$color <- "green"}

  The notation for edges is similar: \code{E(g)} means all edges
  \code{E(g)$weight} is the \code{weight} attribute for all edges, etc.

  See also the manual page for \code{iterators} about how to create
  various vertex and edge sets.
}
\value{
  \code{get.graph.attribute}, \code{get.vertex.attribute} and
  \code{get.edge.attribute} return an R object, or a list of R objects
  if attributes of more vertices/edges are requested.

  \code{set.graph.attribute}, \code{set.vertex.attribute},
  \code{set.edge.attribute}, and also \code{remove.graph.attribute},
  \code{remove.vertex.attribute} and \code{remove.edge.attribute} return
  a new graph object with the updates/removes performed.

  \code{list.graph.attributes}, \code{list.vertex.attributes} and
  \code{list.edge.attributes} return a character vector, the names of
  the attributes present.
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{print.igraph}} can also print attributes}
\examples{
g <- graph.ring(10)
g <- set.graph.attribute(g, "name", "RING")
# It is the same as
g$name <- "RING"
g$name

g <- set.vertex.attribute(g, "color", value=c("red", "green"))
get.vertex.attribute(g, "color")

g <- set.edge.attribute(g, "weight", value=runif(ecount(g)))
get.edge.attribute(g, "weight")

# The following notation is more convenient

g <- graph.star(10)

V(g)$color <- c("red", "green")
V(g)$color

E(g)$weight <- runif(ecount(g))
E(g)$weight

print(g, g=TRUE, v=TRUE, e=TRUE)
}
\keyword{graphs}
