\name{hsbm.game}
\alias{hsbm.game}
\concept{Stochastic block model}
\concept{Random graph model}
\title{Sample hierarchical stochastic block model}
\description{Sampling from a hierarchical stochastic block model of networks}
\usage{
hsbm.game (n, m, rho, C, p)
}
\arguments{
  \item{n}{Integer scalar, the number of vertices.}
  \item{m}{Integer scalar, the number of vertices per block. \code{n /
      m} must be integer.}
  \item{rho}{Numeric vector, the fraction of vretices per cluster,
    within a block. Must sum up to 1, and \code{rho * m} must be integer
    for all elements of rho.}
  \item{C}{A square, symmetric numeric matrix, the Bernoulli rates for
    the clusters within a block. Its size must mach the size of the
    \code{rho} vector.}
  \item{p}{Numeric scalar, the Bernoulli rate of connections between
    vertices that in different blocks.}
}
\details{
  TODO
}
\value{
  An igraph graph.
}
\references{TODO}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{sbm.game}}}
\examples{
## Ten blocks with three clusters each
C <- matrix(c(1  , 3/4,   0,
              3/4,   0, 3/4,
              0  , 3/4, 3/4), nrow=3)
g <- hsbm.game(100, 10, rho=c(3,3,4)/10, C=C, p=1/20)
g
if (require(Matrix)) {
  image(g[])
}
}
\keyword{graphs}
