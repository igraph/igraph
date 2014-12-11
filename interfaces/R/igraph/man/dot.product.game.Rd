\name{sample_dot_product}
\alias{sample_dot_product}
\alias{dot_product}
\title{Generate random graphs according to the random dot product graph
  model}
\description{In this model, each vertex is represented by a latent
  position vector. Probability of an edge between two vertices are given
  by the dot product of their latent position vectors.}
\usage{
sample_dot_product (vecs, directed=FALSE)
dot_product (...)
}
\arguments{
  \item{vecs}{A numeric matrix in which each latent position vector is
    a column. }
  \item{directed}{A logical scalar, TRUE if the generated graph
    should be directed.}
  \item{\dots}{Passed to \code{sample_dot_product}.}
}
\details{
  The dot product of the latent position vectors should be in the [0,1]
  interval, otherwise a warning is given. For negative dot products, no
  edges are added; dot products that are larger than one always add an
  edge.
}
\value{
  An igraph graph object which is the generated random dot product graph.
}
\references{Christine Leigh Myers Nickel: Random dot product graphs, a
  model for social networks. Dissertation, Johns Hopkins University,
  Maryland, USA, 2006.
}
\author{ Gabor Csardi \email{csardi.gabor@gmail.com} }
\seealso{\code{\link{sample_dirichlet}},
  \code{\link{sample_sphere_surface}} and
  \code{\link{sample_sphere_volume}} for sampling position vectors.
}
\examples{
## A randomly generated  graph
lpvs <- matrix(rnorm(200), 20, 10)
lpvs <- apply(lpvs, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
g <- sample_dot_product(lpvs)
g

## Sample latent vectors from the surface of the unit sphere
lpvs2 <- sample_sphere_surface(dim=5, n=20)
g2 <- sample_dot_product(lpvs2)
g2
}
\keyword{graphs}
