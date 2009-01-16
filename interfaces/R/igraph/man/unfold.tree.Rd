\name{unfold.tree}
\alias{unfold.tree}
\concept{Tree}
\concept{Forest}
\concept{Breadth-first search}
\title{Convert a general graph into a forest}
\description{
Perform a breadth-first search on a graph and convert it into a tree or
forest by replicating vertices that were found more than once.
}
\usage{
unfold.tree(graph, mode = c("all", "out", "in", "total"), roots) 
}
\arguments{
  \item{graph}{The input graph, it can be either directed or
    undirected.}
  \item{mode}{Character string, defined the types of the paths used for
    the breadth-first search. \dQuote{out} follows the outgoing,
    \dQuote{in} the incoming edges, \dQuote{all} and \dQuote{total} both
    of them. This argument is ignored for undirected graphs.}
  \item{roots}{A vector giving the vertices from which the breadth-first
    search is performed. Typically it contains one vertex per component.}
}
\details{
  A forest is a graph, whose components are trees.
  
  The \code{roots} vector can be calculated by simply doing a
  topological sort in all components of the graph, see the examples
  below.
}
\value{
  A list with two components:
  \item{tree}{The result, an \code{igraph} object, a tree or a forest.}
  \item{vertex_index}{A numeric vector, it gives a mapping from the
    vertices of the new graph to the vertices of the old graph.}
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% \seealso{}
\examples{
g <- graph.tree(10) %du% graph.tree(10)
V(g)$id <- seq_len(vcount(g))-1
roots <- sapply(decompose.graph(g), function(x) {
            V(x)$id[ topological.sort(x)[1]+1 ] })
tree <- unfold.tree(g, roots=roots)
}
\keyword{graphs}

