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
  \item{p}{a numeric scalar}
  \item{permutation}{TODO}  
  \item{n}{TODO}  
  \item{directed}{TODO}    
  \item{A}{TODO}  
  \item{P}{TODO}  
}
\details{
  TODO
}
\value{
  TODO 
  
}
\references{TODO}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{erdos.renyi.game}}}
\examples{
 g <- erdos.renyi.game(10, .1)
 g2 <- correlated.game(g, corr=1, p=g$p, perm=NULL)
 g3 <- correlated.pair.game(n=10,corr=.8, p=.5,directed=FALSE)
 g3$graph1
 g3$graph2
}
\keyword{graphs}
