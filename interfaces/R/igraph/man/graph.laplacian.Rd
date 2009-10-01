\name{graph.laplacian}
\alias{graph.laplacian}
\concept{Graph Laplacian}
\title{Graph Laplacian}
\description{The Laplacian of a graph.}
\usage{
graph.laplacian(graph, normalized=FALSE)
}
\arguments{
  \item{graph}{The input graph.}
  \item{normalized}{Whether to calculate the normalized Laplacian. See
    definitions below.}
}
\details{
  The Laplacian Matrix of a graph is a symmetric matrix having
  the same number of rows and columns as the number of vertices in the
  graph and element (i,j) is d[i], the degree of vertex i if if i==j,
  -1 if i!=j and there is an edge between vertices i and j and 0
  otherwise.

  A normalized version of the Laplacian Matrix is similar: element (i,j)
  is 1 if i==j, -1/sqrt(d[i] d[j]) if i!=j and there is an edge between
  vertices i and j and 0 otherwise.
}
\value{A square matrix with as many rows as the number of vertices in
  the input graph.}
% \references{}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
g <- graph.ring(10)
graph.laplacian(g)
graph.laplacian(g, norm=TRUE)
}
\keyword{graphs}
