\name{k.regular.game}
\alias{k.regular.game}
\concept{Random graph model}
\concept{Regular graph}
\title{Create a random regular graph}
\description{Generate a random graph where each vertex has the same
  degree.}
\usage{
k.regular.game (no.of.nodes, k, directed = FALSE, multiple = FALSE) 
}
\arguments{
  \item{no.of.nodes}{Integer scalar, the number of vertices in the
    generated graph.}
  \item{k}{Integer scalar, the degree of each vertex in the graph, or
    the out-degree and in-degree in a directed graph.}
  \item{directed}{Logical scalar, whether to create a directed graph.}
  \item{multiple}{Logical scalar, whether multiple edges are allowed.}
}
\details{
  This game generates a directed or undirected random graph where the
  degrees of vertices are equal to a predefined constant k. For undirected
  graphs, at least one of k and the number of vertices must be even.

  The game simply uses \code{\link{degree.dequence.game}} with
  appropriately constructed degree sequences.
}
\value{
  An igraph graph.
}
% \references{}
\author{
  Tamas Nepusz \email{ntamas@gmail.com}
}
\seealso{
  \code{\link{degree.sequence.game}} for a generator with prescribed
  degree sequence.
}
\examples{
## A simple ring
ring <- k.regular.game(10, 2)
plot(ring)

## k-regular graphs on 10 vertices, with k=1:9
k10 <- lapply(1:9, k.regular.game, no.of.nodes=10)

layout(matrix(1:9, nrow=3, byrow=TRUE))
sapply(k10, plot, vertex.label=NA)
}
\keyword{graphs}


