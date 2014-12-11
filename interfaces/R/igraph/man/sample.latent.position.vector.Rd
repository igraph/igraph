\name{sample.latent.position.vector}
\alias{sample.latent.position.vector}
\alias{sample_dirichlet}
\alias{sample_sphere_surface}
\alias{sample_sphere_volume}
\concept{Random vector}
\concept{Random graph model}
\title{Sample finite-dimensional vectors}
\description{Sample finite-dimensional vectors to use as latent position
  vectors in random dot product graphs}
\usage{
sample_dirichlet (n, alpha) 
sample_sphere_surface (dim, n = 1, radius = 1, positive = TRUE) 
sample_sphere_volume  (dim, n = 1, radius = 1, positive = TRUE) 
}
\arguments{
  \item{n}{Integer scalar, the sample size.}
  \item{alpha}{Numeric vector, the vector of \eqn{\alpha}{alpha}
    parameter for the Dirichlet distribution.} 
  \item{dim}{Integer scalar, the dimension of the random vectors.}
  \item{radius}{Numeric scalar, the radius of the sphere to sample.}
  \item{positive}{Logical scalar, whether to sample from the positive
    orthant of the sphere.}
}
\details{
  \code{sample_dirichlet} generates samples from the Dirichlet
  distribution with given \eqn{\alpha}{alpha} parameter. The sample is
  drawn from \code{length(alpha)-1}-simplex.
  
  \code{sample_sphere_surface} generates uniform samples from
  \eqn{S^{dim-1}} (the \code{(dim-1)}-sphere) with radius \code{radius},
  i.e. the Euclidean norm of the samples equal \code{radius}. 
  
  \code{sample_sphere_volume} generates uniform samples from
  \eqn{S^{dim-1}} (the \code{(dim-1)}-sphere) i.e. the Euclidean norm of
  the samples is smaller or equal to  \code{radius}.
 } 
\value{
  A \code{dim} (length of the \code{alpha} vector for
  \code{sample_dirichlet}) times \code{n} matrix, whose columns are the
  sample vectors.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{sample_dot_product}}}
\examples{
lpvs.dir    <- sample_dirichlet(n=20, alpha=rep(1, 10))
RDP.graph.2 <- sample_dot_product(lpvs.dir)
colSums(lpvs.dir)

lpvs.sph    <- sample_sphere_surface(dim=10, n=20, radius=1)
RDP.graph.3 <- sample_dot_product(lpvs.sph)
vec.norm    <- apply(lpvs.sph, 2, function(x) { sum(x^2) })
vec.norm

lpvs.sph.vol <- sample_sphere_volume(dim=10, n=20, radius=1)
RDP.graph.4  <- sample_dot_product(lpvs.sph.vol)
vec.norm     <- apply(lpvs.sph.vol, 2, function(x) { sum(x^2) })
vec.norm
}
\keyword{graphs}
