\name{sample.latent.position.vector}
\alias{sample.latent.position.vector}
\alias{sample.dirichlet}
\alias{sample.sphere.surface}
\alias{sample.sphere.volume}
\concept{Random vector}
\concept{Random graph model}
\title{Sample finite-dimensional vectors  }
\description{Sample finite-dimensional vectors to use as latent position vectors in random dot product graphs}
\usage{
sample.dirichlet (n, alpha) 
sample.sphere.surface (dim, n = 1, radius = 1, positive = TRUE) 
sample.sphere.volume  (dim, n, radius = 1, positive = TRUE) 
}
\arguments{
  \item{n}{integer, the sample size}
  \item{alpha}{numeric vector, the vector of \eqn{\alpha}{alpha} parameter for the Dirichlet distribution }
  \item{dim}{integer, the dimension of the random vectors}
  \item{radius}{numeric, the radius of the sampled sphere. Determines 
  the magnitude of the random vectors}
  \item{positive}{logical , TODO}
  
}
\details{
  TODO
}
\value{
  TODO  
}
\references{TODO}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{sbm.game}}}
\examples{

lpvs.dir <- sample.dirichlet(n=20,rep(1,10))
RDP.graph.2<- dot.product.game(lpvs.dir)
colSums(lpvs.dir)

lpvs.sph <- sample.sphere.surface(dim=10,n=20,radius=1)
RDP.graph.3<- dot.product.game(lpvs.sph)
vec.norm <-apply(lpvs.sph,2,function(x) {sum(x^2)})
vec.norm

lpvs.sph.vol <- sample.sphere.volume(dim=10,n=20,radius=1)
RDP.graph.4<- dot.product.game(lpvs.sph.vol)
vec.norm <- apply(lpvs.sph.vol,2,function(x) {sum(x^2)})
vec.norm
}
\keyword{graphs}
