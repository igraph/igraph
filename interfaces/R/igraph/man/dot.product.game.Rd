\name{dot.product.game}
\alias{dot.product.game}
\title{Generate random graphs according to the random dot product graph
  model}
\description{In this model, each vertex is represented by a latent
  position vector. Probability of an edge between two vertices are given
  by the dot product of their latent position vectors.}
\usage{
dot.product.game (vecs, directed=FALSE)
}
\arguments{
  \item{vecs}{A numeric matrix in which each  latent position vector is
    a column. }
  \item{directed}{A logical variable, TRUE if  the generated graph
    should be directed.}
}
\details{
  The dot product of the latent position vectors should be non-negative.
}
\value{
  An igraph graph which is the randomly generated dot product graph.
}
\references{Sussman, D.L., Tang, M., Fishkind, D.E., Priebe, C.E.
  A Consistent Adjacency Spectral Embedding for Stochastic Blockmodel
  Graphs,  \emph{Journal of the American Statistical Association},
  Vol. 107(499), 2012
}
\author{ Gabor Csardi \email{csardi.gabor@gmail.com} }
\seealso{\code{\link{adjacency.spectral.embedding}},
        \code{\link{correlated.game}}, \code{\link{correlated.game}}}
\examples{
## A randomly generated  graph
lpvs <- matrix(rnorm(200),20,10)
lpvs <- apply(lpvs, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
RDP.graph <- dot.product.game(lpvs)
RDP.graph
}
\keyword{graphs}