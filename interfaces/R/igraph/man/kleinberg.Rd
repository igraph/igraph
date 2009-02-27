\name{kleinberg}
\alias{authority.score}
\alias{hub.score}
\concept{Hub and authority score}
\title{Kleinberg's centrality scores.}
\description{Kleinberg's hub and authority scores.}
\usage{
authority.score (graph, scale = TRUE, options = igraph.arpack.default) 
hub.score (graph, scale = TRUE, options = igraph.arpack.default) 
}
\arguments{
  \item{graph}{The input graph.}
  \item{scale}{Logical scalar, whether to scale the result to have a
    maximum score of one. If no scaling is used then the result vector
    has unit length in the Euclidean norm.}
  \item{options}{A named list, to override some ARPACK options. See
    \code{\link{arpack}} for details.}
}
\details{
  The authority scores of the vertices are defined as the principal
  eigenvector of \eqn{A^T A}{t(A)*A}, where \eqn{A} is the adjacency
  matrix of the graph.

  The hub scores of the vertices are defined as the principal
  eigenvector of \eqn{A A^T}{A*t(A)}, where \eqn{A} is the adjacency
  matrix of the graph.

  Obviously, for undirected matrices the adjacency matrix is symmetric
  and the two scores are the same.
}
\value{
  A named list with members:
  \item{vector}{The authority/hub scores of the vertices.}
  \item{value}{The corresponding eigenvalue of the calculated
    principal eigenvector.}
  \item{options}{Some information about the ARPACK computation, it has
    the same members as the \code{options} member returned by
    \code{\link{arpack}}, see that for documentation.}
}
\references{
  J. Kleinberg. Authoritative sources in a hyperlinked
  environment. \emph{Proc. 9th ACM-SIAM Symposium on Discrete Algorithms},
  1998. Extended version in \emph{Journal of the ACM} 46(1999). Also
  appears as IBM Research Report RJ 10076, May 1997.
}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu} }
\seealso{\code{\link{evcent}} for eigenvector centrality,
  \code{\link{page.rank}} for the Page Rank scores. \code{\link{arpack}}
  for the underlining machinery of the computation.
}
\examples{
## An in-star
g <- graph.star(10)
hub.score(g)$vector
authority.score(g)$vector

## A ring
g2 <- graph.ring(10)
hub.score(g2)$vector
authority.score(g2)$vector
}
\keyword{graphs}
