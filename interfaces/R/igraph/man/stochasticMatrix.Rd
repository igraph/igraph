\name{get.stochastic}
\alias{get.stochastic}
\title{Stochastic matrix of a graph}
\description{Retrieves the stochastic matrix of a graph of class
  \code{igraph}.
}
\usage{
get.stochastic(graph, column.wise = FALSE, sparse =
      getIgraphOpt("sparsematrices"))
}
\arguments{
  \item{graph}{The input graph. Must be of class \code{igraph}.}
  \item{column.wise}{If \code{FALSE}, then the rows of the stochastic
    matrix sum up to one; otherwise it is the columns.}
  \item{sparse}{Logical scalar, whether to return a sparse matrix. The
    \code{Matrix} package is needed for sparse matrices.}
}
\details{
  Let \eqn{M} be an \eqn{n \times n}{n x n} adjacency matrix with real
  non-negative entries. Let us define
  \eqn{D = \textrm{diag}(\sum_{i}M_{1i}, \dots, \sum_{i}M_{ni})}{%
    D=diag( sum(M[1,i], i), ..., sum(M[n,i], i) )}
  
  The (row) stochastic matrix is defined as
  \deqn{W = D^{-1}M,}{W = inv(D) M,}
  where it is assumed that \eqn{D} is non-singular.
  Column stochastic matrices are defined in a symmetric way.
}

\value{
 A regular \R matrix or a matrix of class \code{Matrix} if a
 \code{sparse} argument was \code{TRUE}.
}

\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{get.adjacency}}}
\examples{
## g is a large sparse graph
g <- barabasi.game(n = 10^5, power = 2, directed = FALSE)
W <- get.stochastic(g, sparse=TRUE)

## a dense matrix here would probably not fit in the memory
class(W)

## may not be exactly 1, due to numerical errors
max(abs(rowSums(W))-1)
}
\keyword{graphs}
