\name{graph-operators}
\alias{graph.union}
\alias{graph.disjoint.union}
\alias{graph.intersection}
\alias{graph.compose}
\alias{graph.difference}
\alias{graph.complementer}
\alias{\%c\%}
\alias{\%du\%}
\alias{\%m\%}
\alias{\%s\%}
\alias{\%u\%}
\concept{Graph operators}
\title{Graph operators}
\description{Graph operators handle graphs in terms of set theory.}
\usage{
graph.union(\dots)
graph.disjoint.union(\dots)
graph.intersection(\dots)
graph.compose(g1, g2)
graph.difference(big, small)
graph.complementer(graph, loops=FALSE)
x \%c\% y
x \%du\% y
x \%m\% y
x \%s\% y
x \%u\% y
}
\arguments{
  \item{\dots}{Graph objects or lists of graph objects.}
  \item{g1,g2,big,small,graph,x,y}{Graph objects.}
  \item{loops}{Logical constant, whether to generate loop edges.}
}
\details{
  A graph is homogenous binary relation over the set 0, \dots, |V|-1,
  |V| is the number of vertices in the graph. A homogenous binary
  relation is a set of ordered (directed graphs) or unordered
  (undirected graphs) pairs taken from 0, \dots, |V|-1. The functions
  documented here handle graphs as relations.

  \code{graph.union} creates the union of two or more graphs. Ie. only
  edges which are included in at least one graph will be part of the new
  graph. This function can be also used via the \%u\% operator.

  \code{graph.disjoint.union} creates a union of two or more disjoint
  graphs. Thus first the vertices in the second, third, etc. graphs are
  relabeled to have completely disjoint graphs. Then a simple union is
  created. This function can also be used via the \%du\% operator.

  \code{graph.intersection} creates the intersection of two or more
  graphs: only edges present in all graphs will be included.
  The corresponding operator is \%s\%.

  \code{graph.difference} creates the difference of two graphs. Only
  edges present in the first graph but not in the second will be be
  included in the new graph. The corresponding operator is \%m\%.

  \code{graph.complementer} creates the complementer of a graph. Only
  edges which are \emph{not} present in the original graph will be
  included in the new graph.

  \code{graph.compose} creates the composition of two graphs. The new
  graph will contain an (a,b) edge only if there is a vertex c, such
  that edge (a,c) is included in the first graph and (c,b) is included
  in the second graph. The corresponding operator is \%c\%.

  \code{graph.complementer} keeps graph and vertex attriubutes, edge
  attributes are lost. \code{graph.difference} keeps all attributes
  (graph, vertex and edge) of the first graph.
  The other functions do not handle vertex and edge
  attributes, the new graph will have no attributes at all.
}
\value{
A new graph object.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% \seealso{}
\examples{
g1 <- graph.ring(10)
g2 <- graph.star(10, mode="undirected")
graph.union(g1, g2)
graph.disjoint.union(g1, g2)
graph.intersection(g1, g2)
graph.difference(g1, g2)
graph.complementer(g2)
graph.compose(g1, g2)

## graph complementer keeps graph and vertex attributes
g2$layout <- layout.circle
V(g2)$name <- letters[1:vcount(g2)]
cg2 <- graph.complementer(g2)
cg2$layout
V(cg2)$name

## graph difference keeps the attributes of the first graph
g1$name <- "Ring"
V(g1)$name <- LETTERS[1:vcount(g1)]
E(g1)$weight <- round(runif(ecount(g1))*10)
dg <- graph.difference(g1, g2)
print(dg, g=TRUE, v=TRUE, e=TRUE)
}
\keyword{graphs}
