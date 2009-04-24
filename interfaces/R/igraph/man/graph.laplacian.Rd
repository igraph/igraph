\name{graph.laplacian}
\alias{graph.laplacian}
\concept{Graph Laplacian}
\title{Graph Laplacian}
\description{The Laplacian of a graph.}
\usage{
graph.laplacian(graph, normalized=FALSE, weights=NULL)
}
\arguments{
  \item{graph}{The input graph.}
  \item{normalized}{Whether to calculate the normalized Laplacian. See
    definitions below.}
  \item{weights}{An optional vector giving edge weights for weighted
    Laplacian matrix. If this is \code{NULL} and the graph has an edge
    attribute called \code{weight}, then it will be used
    automatically. Set this to \code{NA} if you want the unweighted
    Laplacian on a graph that has a \code{weight} edge attribute.}
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

  The weighted version of the Laplacian simply works with the weighted
  degree instead of the plain degree. I.e. (i,j) is d[i], the weighted
  degree of vertex i if if i==j, -w if i!=j and there is an edge between
  vertices i and j with weight w, and 0 otherwise. The weighted degree
  of a vertex is the sum of the weights of its adjacent edges.
}
\value{A numeric matrix.}
% \references{}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
g <- graph.ring(10)
graph.laplacian(g)
graph.laplacian(g, norm=TRUE)
}
\keyword{graphs}
