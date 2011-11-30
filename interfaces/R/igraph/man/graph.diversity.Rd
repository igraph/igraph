\name{graph.diversity}
\alias{graph.diversity}
\concept{Entropy}
\title{Graph diversity}
\description{
  Calculates a measure of diversity for all vertices.
}
\usage{
graph.diversity(graph, weights = NULL)
}
\arguments{
  \item{graph}{The input graph. Edge directions are ignored.}
  \item{weights}{\code{NULL}, or the vector of edge weights to use for
    the computation. If \code{NULL}, then the \sQuote{weight} attibute
    is used. Note that this measure is not defined for unweighted
    graphs.}
}
\details{
  The diversity of a vertex is defined as the (scaled) Shannon entropy
  of the weights of its incident edges:
  \deqn{D(i)=\frac{H(i)}{\log k_i}}{D(i)=H(i)/log(k[i])}
  and
  \deqn{H(i)=-\sum_{j=1}^{k_i} p_{ij}\log p_{ij},}{%
    H(i) = -sum(p[i,j] log(p[i,j]), j=1..k[i]),}
  where
  \deqn{p_{ij}=\frac{w_{ij}}{\sum_{l=1}^k_i}V_{il},}{%
    p[i,j] = sum(w[i,l], l=1..k[i]),}
  and \eqn{k_i}{k[i]} is the (total) degree of vertex \eqn{i},
  \eqn{w_{ij}}{w[i,j]} is the weight of the edge(s) between vertices
  \eqn{i} and \eqn{j}.

  For vertices with degree less than two the function returns
  \code{NaN}.  
}
\value{
  A numeric vector, its length is the number of vertices.
}
\references{
  Nathan Eagle, Michael Macy and Rob Claxton: Network Diversity and
  Economic Development, \emph{Science} \bold{328}, 1029--1031, 2010.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% \seealso{}
\examples{
g1 <- erdos.renyi.game(20, 2/20)
g2 <- erdos.renyi.game(20, 2/20)
g3 <- erdos.renyi.game(20, 5/20)
E(g1)$weight <- 1
E(g2)$weight <- runif(ecount(g2))
E(g3)$weight <- runif(ecount(g3))
graph.diversity(g1)
graph.diversity(g2)
graph.diversity(g3)
}
\keyword{graphs}

