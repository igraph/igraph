\name{stochasticMatrix}
\alias{stochasticMatrix}
\title{Stochastic matrix of a graph}
\description{Retrieves the stochastic matrix of a graph of class
  \code{igraph}.
}
\usage{
stochasticMatrix(graph, norm = c("row", "col"), ...)
}
\arguments{
  \item{graph}{The input graph. Must be of class \code{igraph}.}
  \item{norm}{Either \dQuote{row} or \dQuote{col}. If set to
    \dQuote{row} the rows of the Laplacian matrix sum up to zero and the
    rows of the stochastic matrix sum up to one; otherwise it is the
    columns.}
  \item{\dots}{Further arguments to be passed to \code{get.adjacency}.}
}
\details{
  If \code{graph} is a weighted graph attention must be paid to use
  non-negative weights. 

  Let \eqn{M} be an \eqn{n \times n}{n x n} adjacency matrix with real
  non-negative entries. Let us define
  \eqn{D = \textrm{diag}(\sum_{i}M_{1i}, \dots, \sum_{i}M_{ni})}{%
    D=diag( sum(M[1,i], i), ..., sum(M[n,i], i) )}
  
  The (row) stochastic matrix is defined as
  \deqn{W = D^{-1}M,}{W = inv(D) M,}
  where it is assumed that \eqn{D} is non-singular.
  Column stochastic and Laplacian matrices are defined in a symmetric
  way.
}

\value{
 A regular \R matrix or a matrix of class \code{Matrix} if a
 \code{sparse} argument was passed to \code{\link{get.adjacency}}.
}

\author{David Morton de Lachapelle \email{david.morton@epfl.ch},
  \email{david.mortondelachapelle@swissquote.ch}}
\seealso{\code{\link{get.adjacency}}}
\examples{
## X is a large sparse graph
g <- barabasi.game(n = 10^5, power = 2, directed = FALSE)
W <- stochasticMatrix(g, norm = "col")
class(W)
## a dense matrix here would probably not fit in the memory
}
\keyword{graphs}