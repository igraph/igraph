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
\details{TODO}
\value{TODO}
\references{
  M. Zhu, and A. Ghodsi (2006). Automatic dimensionality
  selection from the scree plot via the use of profile
  likelihood. \emph{Computational Statistics and Data Analysis},
  Vol. 51,  918-930.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{adjacency.spectral.embedding}}}
\examples{
lpvs<-matrix(rnorm(200),10,20)
lpvs <- apply(lpvs,2,norm<-function(x) {return (abs(x)/sqrt(sum(x^2)))})
RDP.graph <- dot.product.game(lpvs)
dim.chosen <- dimSelect(adjacency.spectral.embedding(RDP.graph,10)$D)
dim.chosen


lpvs.dir <- sample.dirichlet(n=20,rep(1,10))
RDP.graph.2<- dot.product.game(lpvs.dir)
dim.chosen <- dimSelect(adjacency.spectral.embedding(RDP.graph.2,10)$D)
dim.chosen

lpvs.sph <- sample.sphere.surface(dim=10,n=20,radius=1)
RDP.graph.3<- dot.product.game(lpvs.sph)
dim.chosen <- dimSelect(adjacency.spectral.embedding(RDP.graph.3,10)$D)
dim.chosen

}
\keyword{graphs, embedding, dimensionality selection}