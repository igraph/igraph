\name{permute.vertices}
\alias{permute.vertices}
\concept{Permutation}
\title{Permute the vertices of a graph}
\description{Create a new graph, by permuting vertex ids.}
\usage{
permute.vertices(graph, permutation)
}
\arguments{
  \item{graph}{The input graph, it can directed or undirected.}
  \item{permutation}{A numeric vector giving the permutation to
    apply. The first element is the new id of vertex 0, etc. Every
    number between zero and \code{vcount(graph)-1} must appear exactly
    once.}
}
\details{
  This function creates a new graph from the input graph by permuting
  its vertices according to the specified mapping. Call this function
  with the output of \code{\link{canonical.permutation}} to create
  the canonical form of a graph.  
}
\value{A new graph object.}
%\references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{canonical.permutation}}}
\examples{
# Random permutation of a random graph
g <- random.graph.game(20, 50, type="gnm")
g2 <- permute.vertices(g, sample(vcount(g))-1)
graph.isomorphic(g, g2)
}
\keyword{graphs}
