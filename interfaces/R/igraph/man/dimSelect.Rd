\name{dimSelect}
\alias{dimSelect}
\concept{Dimensionality selection for Spectral Embedding}
\title{Dimensionality selection for Spectral Embedding using profile likelihood}
\description{ }
\usage{
dimSelect(sv)
}
\arguments{
  \item{sv}{}
}
\details{
  }
\value{}
\references{
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{adjacency.spectral.embedding}}}
\examples{
lpvs<-matrix(rnorm(200),20,10)
lpvs <- apply(lpvs,2,norm<-function(x) {return (abs(x)/sqrt(sum(x^2)))})
RDP.graph <- dot.product.game(lpvs)
dim.chosen <- dimSelect(adjacency.spectral.embedding(RDP.graph,9)$D)
dim.chosen
}
\keyword{graphs}
