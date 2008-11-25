\name{graphMatrices}
\alias{adjacencyMatrix}
\alias{laplacianMatrix}
\alias{stochasticMatrix}
\title{Common Graph Matrices}
\description{These functions retrieve the adjacency, the Laplacian and the stochastic matrix of a graph of class \code{igraph}. The last two functions can also convert a regular \R matrix or a sparse matrix (of class \code{Matrix}) into a Laplacian, respectively a stochastic matrix.
}
\usage{
adjacencyMatrix(X, sparse = TRUE, ...)
laplacianMatrix(X, norm = c("row", "col"), sparse = TRUE, 
		normalized = FALSE, ...)
stochasticMatrix(X, norm = c("row", "col"), sparse = TRUE, ...)
}
\arguments{
\item{X}{The input graph or square matrix. Must be of class \code{igraph}, \code{matrix} or \code{Matrix}.}

\item{sparse}{When \code{X} is a graph, the \code{sparse} argument passed on to \code{get.adjacency} of the \pkg{igraph} package. 
Must be \code{TRUE} or \code{FALSE}.}

\item{norm}{Either \dQuote{row} or \dQuote{col}. If set to \dQuote{row} the rows of the Laplacian matrix sum up to zero and the rows 
of the stochastic matrix sum up to one; otherwise it is the columns.}

\item{normalized}{Either \code{TRUE} or \code{FALSE}. If \code{TRUE} returns with the normalized Laplacian (see Details).}

\item{\dots}{Further arguments to be passed to \code{get.adjacency} (only when \code{X} is a graph).}
}
\details{
If \code{X} is a square matrix it must have non-negative entries for the stochastic and the Laplacian matrices to be defined. Likewise, if \code{X} is a weighted graph attention must be paid to use non-negative weights. 

The three functions call \code{get.adjacency} of the \pkg{igraph} package, which in turn requires the package \pkg{Matrix} if \code{sparse} equals \code{TRUE}. \pkg{igraph} and \pkg{Matrix} should be installed along with the  \pkg{SCG} package. 

Let \eqn{M} be a \eqn{n \times n} real matrix with non-negative entries. The (row) Laplacian matrix \eqn{L} is defined as
	\deqn{L = D - M,}
where \eqn{D = \textrm{diag}(\sum_{i}M_{1i}, \dots, \sum_{i}M_{ni}) } .
The normalized Laplacian matrix, introduced by Fan R. K. Chung in \cite{Spectral Graph Theory, AMS, 1994}, is defined as
	\deqn{L_N = D^{-1/2}LD^{-1/2},}
where \eqn{L} is the (row) Laplacian Matrix with the extra convention that \eqn{D_{ii}^{-1}=0} if \eqn{\sum_{k}M_{ik} = 0}.
The (row) stochastic matrix is defined as
\deqn{W = D^{-1}M,}
where it is assumed that \eqn{D} is non-singular.
Column stochastic and Laplacian matrices are defined in a symmetric way.
}

\value{
 A regular \R matrix or a matrix of class \code{Matrix} if \code{sparse} = \code{TRUE} (default).
}

\author{David Morton de Lachapelle \email{david.morton@epfl.ch}, \email{david.mortondelachapelle@swissquote.ch}}
\seealso{\code{\link{get.adjacency}}}

\examples{
# X is a random matrix
m <- matrix(runif(100),10,10)
rowSums(m)
L <- laplacianMatrix(m, sparse = FALSE)
rowSums(L)

## X is a large sparse graph
g <- barabasi.game(n = 10^5, power = 2, directed = FALSE)
W <- stochasticMatrix(g, norm = "col")
class(W)
## a dense matrix here would probably not fit in the memory
}

\keyword{array}
\keyword{graphs}