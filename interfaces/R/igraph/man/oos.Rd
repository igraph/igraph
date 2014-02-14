\name{oos}
\alias{oos}
\concept{Out-of-sample embedding}
\title{Out-of-sample embedding for adjacency spectral embedding}
\description{Performs out-of-sample embedding given an existing embedding for the in-sample points and connectivity pattern between the in-sample points and the out-of-sample points}
\usage{
oos(A,X) 
}
\arguments{
\item{A}{A rectangular adjacency matrix of size nxm possibly indexed from an igraph graph object}
\item{X}{A nxd matrix of embeddings for the in-sample points in R^d}
}
\value{
A matrix of size mxd for the embeddings of the out-of-sample points.
}

\author{ Minh Tang \email{mtang10@jhu.edu}}

\examples{
###### sample usage  

n <- 200
m <- 100

X <- matrix(runif(2*(n+m), max = 1/sqrt(2)), nrow = 2)
g <- dot.product.game(X)
g.insample <- induced.subgraph(g, 1:n)
Xhat.insample <- adjacency.spectral.embedding(g.insample, 2)$X
Xhat.oos <- oos(g[][1:n,(n+1):(n+m)], Xhat.insample)
}
