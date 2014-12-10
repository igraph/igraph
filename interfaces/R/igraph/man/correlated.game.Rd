\name{correlated.game}
\alias{correlated.game}
\alias{correlated.pair.game}
\concept{Correlated Graph Pair model}
\concept{Random graph model}
\title{Generate a new random graph from a given graph by randomly
  adding/removing edges.}
\description{Sample a new graph by perturbing the adjacency matrix of a
  given graph and shuffling its vertices.}
\usage{
correlated.game (old.graph, corr, p = old.graph$p, permutation = NULL) 
correlated.pair.game (n, corr, p, directed = FALSE, permutation = NULL)
}
\arguments{
  \item{old.graph}{The original graph.}
  \item{corr}{A scalar in the unit interval, the target Pearson
    correlation between the adjacency matrices of the original the
    generated graph (the adjacency matrix being used as a vector).}
  \item{p}{A numeric scalar, the probability of an edge between two
    vertices, it must in the open (0,1) interval.}
  \item{permutation}{A numeric vector, a permutation vector that is
    applied on the vertices of the first graph, to get the second graph.
    If \code{NULL}, the vertices are not permuted.}
  \item{n}{Numeric scalar, the number of vertices for the sampled
    graphs.}
  \item{directed}{Logical scalar, whether to generate directed graphs.}
}  
\details{
  Please see the reference given below.
}
\value{
  \code{correlated.game} generates an unweighted graph  of the same size as 
  \code{old.graph} such that the correlation coefficient between the entries
  of the two adjacency matrices is \code{corr}. 
  Note  each pair of corresponding matrix entries is a pair of  
  correlated Bernoulli random variables.
  
  \code{correlated.pair.game} generates a list of two igraph objects, 
  named \code{graph1} and \code{graph2}, which are two graphs 
  whose adjacency matrix entries are correlated with \code{corr}.
  
}
\references{Lyzinski, V., Fishkind, D. E., Priebe, C. E. (2013).
  Seeded graph matching for correlated Erdos-Renyi graphs.
  \url{http://arxiv.org/abs/1304.7844}
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{erdos.renyi.game}}}
\examples{
g <- erdos.renyi.game(1000, .1)
g2 <- correlated.game(g, corr=0.5)
cor(as.vector(g[]), as.vector(g2[]))
g
g2

gg <- correlated.pair.game(n=10, corr=.8, p=.5, directed=FALSE)
gg
}
\keyword{graphs,random graphs}
