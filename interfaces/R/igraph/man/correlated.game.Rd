\name{correlated.game}
\alias{correlated.game}
\alias{correlated.pair.game}
\alias{adjcorr}
\concept{Correlated Graph Pair model}
\concept{Random graph model}
\title{Generate a new random graph from a given graph by randomly adding/removing edges }
\description{Sample a new graph by perturbing the adjacency matrix of a given graph and shuffling its vertices}
\usage{
correlated.game (old.graph, corr, p = old.graph$p, permutation = NULL) 
correlated.pair.game (n, corr, p, directed = FALSE, permutation = NULL)
adjcorr(A, P, corr, permutation)
}
\arguments{
  \item{old.graph}{The original graph}
  \item{corr}{a scalar}
  \item{p}{a numeric, the probability of an edge between two vertices}
  \item{permutation}{if NULL, old.graph and the new graph (graph1 and graph2) }  
  \item{n}{numeric, the number of vertices for the sampled graph}  
  \item{directed}{logical, if TRUE, the adjacency matrices 
  are not necessarily symmetric matrices and each entry of the matrices
  is independent of the other entries. 
  If FALSE, the adjacency matrices are constrained to be symmetric 
  and the same Bernoulli random variables determine both upper 
  and lower-diagonal entries.}    
  \item{A}{the adjacency matrix of the first graph}  
  \item{P}{TODO}  
}  
\details{
  TODO
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
Online. \url{http://arxiv.org/abs/1304.7844}}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{erdos.renyi.game}}}
\examples{
 g <- erdos.renyi.game(10, .1)
 g2 <- correlated.game(g, corr=1, p=g$p, perm=NULL)
 g3 <- correlated.pair.game(n=10, corr=.8, p=.5, directed=FALSE)
 g3$graph1
 g3$graph2
}
\keyword{graphs,random graphs}
