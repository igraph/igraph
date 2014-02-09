\name{dimSelect}
\alias{dimSelect}
\concept{Dimensionality selection for Spectral Embedding}
\title{Dimensionality selection for Spectral Embedding using
profile  likelihood.}
\description{For a graph to be embedded using Spectral Embedding,
this function returns an embedding dimension for spectral embedding.
Its only argument is the list of singular values of the adjacency matrix.
The embedding dimension choice is based on the \cite{Zhu & Ghodsi} paper
on automatic dimensionality selection.}
\usage{
dimSelect(sv)
}
\arguments{
  \item{sv}{a numeric vector, the singular values of the adjacency matrix
  of the graph}
}
\details{The input of the function is a numeric vector which contains the measure of \sQuote{importance} for each dimension. For spectral embedding, these are the singular values of the adjacency matrix. The singular values are assumed to be generated from a Gaussian mixture distribution with two components that have different means and same variance. The  dimensionality \eqn{d} is chosen to maximize the likelihood when the $d$ largest singular values are assigned to one component of the mixture and the rest of the singular values assigned to the other component.}
\value{a numeric scalar, the estimate of \eqn{d} }
\references{
  M. Zhu, and A. Ghodsi (2006). Automatic dimensionality
  selection from the scree plot via the use of profile
  likelihood. \emph{Computational Statistics and Data Analysis},
  Vol. 51,  918-930.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{adjacency.spectral.embedding}}}
\examples{

# Generate the two groups of singular values with 
# Gaussian mixture of two components that have different means
sing.vals  <- c( rnorm (10, mean=1, sd= 1), rnorm(10, mean=3, sd=1) )
sing.vals  <- sort(sing.vals)
dim.chosen <- dimSelect(sing.vals)
dim.chosen

# Sample random vectors with multivariate  normal distribution
# and normalize to unit length
lpvs <- matrix(rnorm(200), 10, 20)
lpvs <- apply(lpvs, 2, norm<-function(x) {return (abs(x)/sqrt(sum(x^2)))})
RDP.graph  <- dot.product.game(lpvs)
dim.chosen <- dimSelect( adjacency.spectral.embedding(RDP.graph, 10)$D )
dim.chosen

# Sample random vectors with the Dirichlet distribution
lpvs.dir    <- sample.dirichlet(n=20, rep(1, 10))
RDP.graph.2 <- dot.product.game(lpvs.dir)
dim.chosen  <- dimSelect( adjacency.spectral.embedding(RDP.graph.2, 10)$D )
dim.chosen 


# Sample random vectors from hypersphere with radius 1.
lpvs.sph    <- sample.sphere.surface(dim=10, n=20, radius=1)
RDP.graph.3 <- dot.product.game(lpvs.sph)
dim.chosen  <- dimSelect( adjacency.spectral.embedding(RDP.graph.3, 10)$D )
dim.chosen

}
\keyword{graphs, embedding, dimensionality selection}