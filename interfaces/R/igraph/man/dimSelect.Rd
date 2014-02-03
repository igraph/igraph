\name{dimSelect}
\alias{dimSelect}
\concept{Dimensionality selection for Spectral Embedding}
\title{Dimensionality selection for Spectral Embedding using profile
  likelihood.}
\description{TODO}
\usage{
dimSelect(sv)
}
\arguments{
  \item{sv}{TODO}
}
\details{TODO}
\value{TODO}
\references{
  M. Zhu, and A. Ghodsi (2006). Automatic dimensionality
  selection from the scree plot via the use of profile
  likelihood. \emph{Computational Statistics & Data Analysis}, Vol. 51,
  918-930.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{adjacency.spectral.embedding}}}
\examples{
lpvs <- matrix(rnorm(200), 20, 10)
lpvs <- apply(lpvs, 2 , function(x) { return (abs(x)/sqrt(sum(x^2))) })
RDP.graph <- dot.product.game(lpvs)
dimSelect(adjacency.spectral.embedding(RDP.graph, 9)$D)
}
\keyword{graphs}
