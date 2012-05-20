\name{subgraph.centrality}
\alias{subgraph.centrality}
\concept{Subgraph centrality}
\title{Find subgraph centrality scores of network positions}
\description{Subgraph centrality of a vertex measures the number of
  subgraphs a vertex participates in, weighting them according to their
  size.}
\usage{
subgraph.centrality (graph, diag=FALSE)
}
\arguments{
  \item{graph}{The input graph, it should be undirected, but the
    implementation does not check this currently. }
  \item{diag}{Boolean scalar, whether to include the diagonal of the
    adjacency matrix in the analysis. Giving \code{FALSE} here
    effectively eliminates the loops edges from the graph before the
    calculation.}
}
\details{
  The subgraph centrality of a vertex is defined as the number of closed
  loops originating at the vertex, where longer loops are exponentially
  downweighted.

  Currently the calculation is performed by explicitly calculating all
  eigenvalues and eigenvectors of the adjacency matrix of the
  graph. This effectively means that the measure can only be calculated
  for small graphs.
}
\value{
  A numeric vector, the subgraph centrality scores of the vertices.
}
\references{
  Ernesto Estrada, Juan A. Rodriguez-Velazquez: Subgraph centrality in
  Complex Networks. \emph{Physical Review E} 71, 056103 (2005).
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com} based on the Matlab
  code by Ernesto Estrada}
\seealso{\code{\link{evcent}}, \code{\link{page.rank}}}
\examples{
g <- ba.game(100, m=4, dir=FALSE)
sc <- subgraph.centrality(g)
cor(degree(g), sc)
}
\keyword{graphs}

