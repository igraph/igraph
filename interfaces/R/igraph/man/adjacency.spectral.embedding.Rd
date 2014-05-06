\name{adjacency.spectral.embedding}
\alias{adjacency.spectral.embedding}
\title{Spectral Embedding of Adjacency Matrices}
\description{Spectral decomposition of the adjacency matrices of
  graphs.}
\usage{
adjacency.spectral.embedding (graph, no, weights=NULL,
    which=c("la", "lm", "sm"), scaled = TRUE,
    cvec=degree(graph)/(vcount(graph) - 1),
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
  \item{which}{Which eigenvalue to use. \sQuote{lm} means the one with
    the largest magnitude, \sQuote{la} is the one (algebraic) largest,
    and \sQuote{sa} is the (algebraic) smallest eigenvalue. The default
    is \sQuote{la}.}
  \item{scaled}{Logical scalar, if \code{FALSE}, then \eqn{U} and
    \eqn{V} are returned instead of \eqn{X} and \eqn{Y}.}
  \item{cvec}{A numeric vector, its length is the number vertices in the
    graph. This vector is added to the diagonal of the adjacency
    matrix.}
  \item{options}{A named list containing the parameters for the SVD
    computation algorithm in ARPACK. By default, the list of values is
    assigned the values given by \code{\link{igraph.arpack.default}}.}
}
\details{
  This function computes a \code{no}-dimensional Euclidean
  representation of the graph based on its adjacency
  matrix, \eqn{A}. This representation is computed via the singular value
  decomposition of the adjacency matrix, \eqn{A=UDV^T}.In the case,
  where the graph is a random dot product graph generated using latent
  position vectors in \eqn{R^{no}} for each vertex, the embedding will
  provide an estimate of these latent vectors.

  For undirected graphs the latent positions are calculated as
  \eqn{X=U^{no}D^{1/2}}{U[no] sqrt(D[no])}, where \eqn{U^{no}}{U[no]}
  equals to the first \code{no} columns of \eqn{U}, and
  \eqn{D^{1/2}}{sqrt(D[no])} is a diagonal matrix containing the top
  \code{no} singular values on the diagonal.

  For directed graphs the embedding is defined as the pair
  \eqn{X=U^{no}D^{1/2}}{U[no] sqrt(D[no])} and
  \eqn{Y=V^{no}D^{1/2}}{V[no] sqrt(D[no])}. (For undirected graphs
  \eqn{U=V}, so it is enough to keep one of them.)
}
\value{
  A  list containing with entries:
  \item{X}{Estimated latent positions, an \code{n} times \code{no}
    matrix, \code{n} is the number of vertices.}
  \item{Y}{\code{NULL} for undirected graphs, the second half of the
  latent positions for directed graphs, an \code{n} times \code{no}
    matrix, \code{n} is the number of vertices.}
  \item{options}{A named list, information about the underlying ARPACK
    computation. See \code{\link{arpack}} for the details.}
}
\references{Sussman, D.L., Tang, M., Fishkind, D.E., Priebe, C.E.
  A Consistent Adjacency Spectral Embedding for Stochastic Blockmodel
  Graphs,  \emph{Journal of the American Statistical Association},
  Vol. 107(499), 2012
}
\author{ Gabor Csardi \email{csardi.gabor@gmail.com} }
\seealso{\code{\link{dot.product.game}}}
\examples{
## A small graph
lpvs <- matrix(rnorm(200), 20, 10)
lpvs <- apply(lpvs, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
RDP <- dot.product.game(lpvs)
embed <- adjacency.spectral.embedding(RDP, 5)
}
\keyword{graphs}
