\name{graph.constructors}
\alias{graph.empty}
\alias{graph}
\alias{graph.adjacency}
\alias{graph.star}
\alias{graph.lattice}
\alias{graph.ring}
\alias{graph.tree}
\title{Various methods for creating graphs}
\description{}
\synopsis{
graph.lattice(dimvector = NULL, length = NULL, dim = NULL, nei = 1,
              directed = FALSE, mutual = FALSE, circular = FALSE, \dots)
}
\usage{
graph.empty(\dots, type = igraph.par("default.type"))
graph(edges, \dots)
graph.adjacency(adjmatrix, directed = TRUE, \dots)
graph.star(n, mode = "in", center = 1, directed = TRUE, \dots)
graph.lattice(dimvector, nei = 1, directed = FALSE, mutual = FALSE, 
              circular = FALSE, \dots)
graph.lattice(length, dim, nei = 1, directed = FALSE, mutual = FALSE,
              circular = FALSE, \dots)
graph.ring(n, circular = TRUE, directed = FALSE, mutual = FALSE, \dots)
graph.tree(n, children = 2, \dots)
}
\arguments{
  \item{type}{}
  \item{edges}{}
  \item{adjmatrix}{}
  \item{directed}{}
  \item{n}{}
  \item{mode}{}
  \item{dimvector}{}
  \item{nei}{}
  \item{mutual}{}
  \item{circular}{}
  \item{length}{}
  \item{dim}{}
  \item{children}{}
}
\details{}
\value{}
\references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{}
\examples{}
\keyword{graphs}
