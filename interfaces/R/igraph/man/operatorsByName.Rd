\name{graph-operators-by-name}
\alias{graph.intersection.by.name}
\alias{graph.union.by.name}
\alias{graph.difference.by.name}
\concept{Graph operators}
\title{Graph operators based on symbolic vertex names}
\usage{
graph.intersection.by.name(g1, g2, keep.x.vertices = FALSE,
      keep.x.vertex.attributes = FALSE,
      keep.x.edge.attributes = FALSE)
graph.union.by.name(g1, g2)
graph.difference.by.name(g1, g2, keep.x.vertex.attributes = FALSE,
      keep.x.edge.attributes = FALSE)
}
\arguments{
  \item{g1,g2}{The input graphs, they should be both directed or both
    undirected.}
  \item{keep.x.vertices}{Logical scalar, whether to keep all the
    vertices of \code{g1}, even if they do not appear in \code{g2}.}
  \item{keep.x.vertex.attributes}{Logical scalar, whether to keep all
    vertex attributes of the \code{g1}, even if they do not appear in
    \code{g2}.}
  \item{keep.x.edge.attributes}{Logical scalar, whether to keep edge
    attributes of \code{g1}, even if they do not appear in \code{g2}.}
}
\details{
  These functions treat the graphs as sets of ordered (if they are
  directed), or unordered (if they are undirected) pairs of symbolic
  vertex names and perform set operations on them.
}
\value{
  A new graph object.
}
% \references{}
\author{Magnus Torfason}
\seealso{\code{\link{graph.union}}, \code{link{graph.intersection}} and
  \code{link{graph.difference}} for the functions that are based on
  vertex ids.}
\examples{
# Print attributes as well
igraph.par("print.vertex.attributes", TRUE)
igraph.par("print.edge.attributes", TRUE)

# Creating a wheel graph
g1 <- graph.ring(10)
V(g1)$name <- letters[1:vcount(g1)]
g2 <- graph.star(11, mode="undirected")
V(g2)$name <- c("z", letters[1:vcount(g1)])
g <- graph.union.by.name(g1, g2)
if (interactive()) {
  plot(g, layout=layout.kamada.kawai, vertex.label=V(g)$name)
}
g

# Some more examples
g1 <- graph.formula(a-b-c)
V(g1)$v.attr=c(1,2,3)
E(g1)$e.attr=c(5,7)
g2 <- graph.formula(b-c-d)

# Test the functions
graph.intersection.by.name(g1,g2) # Vertices are intersected as well
graph.union.by.name(g1,g2)        # Vertices are unioned as well
graph.difference.by.name(g1,g2)   # Vertices from x (g1) are used

# graph.intersection.by.name() has some extra parameters
graph.intersection.by.name(g1,g2,
    keep.x.vertices          = TRUE) # Keep all x vertices (only intersect edges)

graph.intersection.by.name(g1,g2,
    keep.x.vertices          = FALSE, # Keep all x vertices (only intersect edges)
    keep.x.vertex.attributes = TRUE, # Don't throw away V(g1) attributes
    keep.x.edge.attributes   = TRUE) # Don't throw away E(g1) attributes

# graph.difference.by.name() has some extra parameters
graph.difference.by.name(g1,g2,
    keep.x.vertex.attributes = TRUE, # Don't throw away V(g1) attributes
    keep.x.edge.attributes   = TRUE) # Don't throw away E(g1) attributes
}
\keyword{graphs}
