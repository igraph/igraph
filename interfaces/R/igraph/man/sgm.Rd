\name{sgm}
\alias{sgm}
\alias{seeded.graph.match}
\concept{Seeded Graph Matching}
\title{Match Graphs given a seeding of vertex correspondences}
\description{ Given two adjacency matrices \code{A} and \code{B} of the same size,
match the two graphs with the help of \code{m} seed vertex pairs which correspond 
to the first \code{m} rows (and columns) of the adjacency matrices. }

\usage{
sgm(A, B, m, start, iteration)
}
\arguments{
  \item{A}{a numeric matrix, the adjacency matrix of the first graph}
  \item{B}{a numeric matrix, the adjacency matrix of the second graph}
  \item{m}{The number of seeds. The first \code{m} vertices of both graphs
  are matched.}
  \item{start}{ a numeric matrix, the  permutation matrix estimate is 
  initialized with \code{start} }
  \item{iteration}{The number of iterations for the Frank-Wolfe algorithm}
}
\details{
The approximate graph matching problem is to find a bijection 
between the vertices of two graphs , such that the number of edge disagreements between the corresponding vertex pairs is minimized. For seeded graph matching,  
part of the bijection that consist of known correspondences (the seeds) is known
and the problem task is to complete the bijection by estimating the permutation
matrix that permutes the rows and columns of the adjacency matrix of the second
graph. 

It is assumed that for the two supplied adjacency matrices \code{A} and \code{B},
both  of size \eqn{n\times n}{n*n}, the first \eqn{m} rows(and columns) of 
\code{A} and \code{B} correspond to the same vertices in both graphs. That is,
the \eqn{n \times n}{n*n} permutation matrix that defines the bijection 
is \eqn{I_{m} \bigoplus  P} for a \eqn{(n-m)\times  (n-m)}{(n-m)*(n-m)} 
permutation matrix \eqn{P} and \eqn{m} times \eqn{m} identity matrix \eqn{I_{m}}.
The function  \code{sgm} estimates the permutation matrix \eqn{P} via 
an optimization algorithm  based on the Frank-Wolfe algorithm.

See references for further details.
}
\value{A numeric matrix which is the permutation matrix that determines
the bijection between the graphs of \code{A} and \code{B} 
}
\references{ Vogelstein, J. T., Conroy, J. M., Podrazik, L. J., 
Kratzer, S. G., Harley, E. T., Fishkind, D. E.,Vogelstein, R. J.,
Priebe, C. E. (2011). 
Fast Approximate Quadratic Programming for Large (Brain) Graph Matching.
Online:  \url{http://arxiv.org/abs/1112.5507}

Fishkind, D. E., Adali, S., Priebe, C. E. (2012). Seeded Graph Matching
Online:  \url{http://arxiv.org/abs/1209.0367}
}
\author{Vince Lyzinski \url{http://www.ams.jhu.edu/~lyzinski/}}
\seealso{\code{\link{correlated.game}},\code{\link{correlated.pair.game}}}
\examples{
 #require(Matrix)
 g1 <- erdos.renyi.game(10, .1)
 randperm <- c(1:3, 3+sample(7))
 g2 <- correlated.game(g1, corr=1, p=g1$p, perm=randperm)
 A  <- as.matrix(get.adjacency(g1))
 B  <- as.matrix(get.adjacency(g2))
 P  <-sgm (A, B, m=3, start=diag(rep(1, nrow(A)-3)), 20)
 P
 }
\keyword{graphs}
