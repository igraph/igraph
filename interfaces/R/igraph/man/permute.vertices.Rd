\name{permute}
\alias{permute.vertices}
\alias{permute}
\concept{Permutation}
\title{Permute the vertices of a graph}
\description{Create a new graph, by permuting vertex ids.}
\usage{
permute(graph, permutation)
}
\arguments{
  \item{graph}{The input graph, it can directed or undirected.}
  \item{permutation}{A numeric vector giving the permutation to
    apply. The first element is the new id of vertex 1, etc. Every
    number between one and \code{vcount(graph)} must appear exactly
    once.}
}
\details{
  This function creates a new graph from the input graph by permuting
  its vertices according to the specified mapping. Call this function
  with the output of \code{\link{canonical_permutation}} to create
  the canonical form of a graph.

  \code{permute} keeps all graph, vertex and edge attributes of
  the graph.
}
\value{A new graph object.}
%\references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{canonical_permutation}}}
\examples{
# Random permutation of a random graph
g <- sample_gnm(20, 50)
g2 <- permute(g, sample(vcount(g)))
graph.isomorphic(g, g2)

# Permutation keeps all attributes
g$name <- "Random graph, Gnm, 20, 50"
V(g)$name <- letters[1:vcount(g)]
E(g)$weight <- sample(1:5, ecount(g), replace=TRUE)
g2 <- permute(g, sample(vcount(g)))
graph.isomorphic(g, g2)
g2$name
V(g2)$name
E(g2)$weight
all(sort(E(g2)$weight) == sort(E(g)$weight))
}
\keyword{graphs}
