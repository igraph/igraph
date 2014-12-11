\name{embed_laplacian_matrix}
\alias{embed_laplacian_matrix}
\title{Spectral Embedding of the Laplacian of a Graph}
\description{Spectral decomposition of Laplacian matrices of
  graphs.}
\usage{
embed_laplacian_matrix (graph, no, weights = NULL,
    which = c("lm", "la", "sa"), degmode = c("out", "in", "all", "total"),
    type = c("default", "D-A", "DAD", "I-DAD", "OAP"), scaled = TRUE,
    options = igraph.arpack.default)
}
\arguments{
  \item{graph}{The input graph, directed or undirected.}
  \item{no}{An integer scalar. This value is the embedding dimension of
    the spectral embedding. Should be smaller than the number of
    vertices. The largest \code{no}-dimensional non-zero
    singular values are used for the spectral embedding.}
  \item{weights}{Optional positive weight vector for calculating
    weighted closeness. If the graph has a \code{weight} edge
    attribute, then this is used by default.}
  \item{which}{Which eigenvalues (or singular values, for directed
    graphs) to use. \sQuote{lm} means the ones with
    the largest magnitude, \sQuote{la} is the ones (algebraic) largest,
    and \sQuote{sa} is the (algebraic) smallest eigenvalues. The default
    is \sQuote{lm}. Note that for directed graphs \sQuote{la} and
    \sQuote{lm} are the equivalent, because the singular values are used
    for the ordering.}
  \item{degmode}{TODO}
  \item{type}{The type of the Laplacian to use. Various definitions
    exist for the Laplacian of a graph, and one can choose
    between them with this argument.

    Possible values:
    \code{D-A} means \eqn{D-A} where \eqn{D} is the
    degree matrix and \eqn{A} is the adjacency matrix;
    \code{DAD} means \eqn{D^{1/2}}{D^1/2} times \eqn{A}
    times \eqn{D^{1/2}{D^1/2}}, \eqn{D^{1/2}}{D^1/2} is the inverse of the
    square root of the degree matrix; \code{I-DAD} means
    \eqn{I-D^{1/2}}{I-D^1/2}, where \eqn{I} is the identity matrix.
    \code{OAP} is \eqn{O^{1/2}AP^{1/2}}{O^1/2 A P^1/2}, where
    \eqn{O^{1/2}}{O^1/2} is the inverse of the square root of the
    out-degree matrix and \eqn{P^{1/2}}{P^1/2} is the same for the in-degree
    matrix.

    \code{OAP} is not defined for undireted graphs, and is the only defined
    type for directed graphs.

    The default (i.e. type \code{default}) is to use \code{D-A} for undirected
    graphs and \code{OAP} for directed graphs.}
  \item{scaled}{Logical scalar, if \code{FALSE}, then \eqn{U} and
    \eqn{V} are returned instead of \eqn{X} and \eqn{Y}.}
  \item{options}{A named list containing the parameters for the SVD
    computation algorithm in ARPACK. By default, the list of values is
    assigned the values given by \code{\link{igraph.arpack.default}}.}
}
\details{
  This function computes a \code{no}-dimensional Euclidean
  representation of the graph based on its Laplacian
  matrix, \eqn{L}. This representation is computed via the singular value
  decomposition of the Laplacian matrix.

  They are essentially doing the same as
  \code{\link{embed_adjacency_matrix}}, but work on the Laplacian
  matrix, instead of the adjacency matrix.
}
\value{
  A  list containing with entries:
  \item{X}{Estimated latent positions, an \code{n} times \code{no}
    matrix, \code{n} is the number of vertices.}
  \item{Y}{\code{NULL} for undirected graphs, the second half of the
  latent positions for directed graphs, an \code{n} times \code{no}
    matrix, \code{n} is the number of vertices.}
  \item{D}{The eigenvalues (for undirected graphs) or the singular
    values (for directed graphs) calculated by the algorithm.}
  \item{options}{A named list, information about the underlying ARPACK
    computation. See \code{\link{arpack}} for the details.}
}
\references{Sussman, D.L., Tang, M., Fishkind, D.E., Priebe, C.E.
  A Consistent Adjacency Spectral Embedding for Stochastic Blockmodel
  Graphs,  \emph{Journal of the American Statistical Association},
  Vol. 107(499), 2012
}
\author{ Gabor Csardi \email{csardi.gabor@gmail.com} }
\seealso{\code{\link{embed_adjacency_matrix}},
  \code{\link{sample_dot_product}}}
\examples{
## A small graph
lpvs <- matrix(rnorm(200), 20, 10)
lpvs <- apply(lpvs, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
RDP <- sample_dot_product(lpvs)
embed <- embed_laplacian_matrix(RDP, 5)
}
\keyword{graphs}
