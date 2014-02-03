\name{adjacency.spectral.embedding}
\alias{adjacency.spectral.embedding}
\title{Spectral Embedding of Adjacency Matrices}
\description{Using a spectral decomposition of the adjacency matrices of graphs}
\usage{
adjacency.spectral.embedding (graph, no,cvec=degree(graph)/(vcount(graph) - 1), 
    options = igraph.arpack.default)
}
\arguments{
  \item{graph}{The input graph. It might be directed, but edge
    directions are ignored.}
  \item{no}{An integer scalar. This value is the embedding dimension of the spectral embedding.. Should be smaller than the number of vertices, in which case the largest \code{no}-dimensional non-zero  singular values are used for the spectral embedding.}   
  \item{cvec}{A numeric vector}
  \item{options}{ A list containing the parameters for the SVD computation algorithm in ARPACK. By default, the list of values is assigned the values given by \code{igraph.arpack.default }.}
 
}
\details{
  This function computes a \code{no}-dimensional Euclidean representation of the graph based on its adjacency matrix,\eqn{A}. This representation is computed via the singular value decomposition of the adjacency matrix, \eqn{A=UDV^T}.In the case, where the graph is a random dot product graph generated using latent position vectors in \eqn{R^{no}} for each vertex , the embedding will provide an estimate of  these latent vectors.
}
\value{
  A  list containing  the matrices D(a \code{no} \eqn{\times}  \code{no} matrix), U(\code{vcount(graph)} \eqn{\times}  no ) and V (\code{vcount(graph)} \ \eqn{\times}  no ) respectively.
}
\references{Sussman, D.L., Tang, M., Fishkind, D.E., \& Priebe, C.E.
A Consistent Adjacency Spectral Embedding for Stochastic Blockmodel Graphs,  \emph{Journal of the American Statistical Association}, Vol. 107(499), 2012
}
\author{ Gabor Csardi \email{csardi.gabor@gmail.com} }
\seealso{\code{\link{dot.product.game}}}
\examples{
## A small graph
lpvs<-matrix(rnorm(200),20,10)
lpvs <- apply(lpvs,2,norm<-function(x) {return (abs(x)/sqrt(sum(x^2)))})
RDP.graph <- dot.product.game(lpvs)
graph.as.embed.config<-adjacency.spectral.embedding(RDP.graph,5)

}
\keyword{graphs,embedding}
