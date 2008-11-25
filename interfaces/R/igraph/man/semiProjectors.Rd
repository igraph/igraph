\name{semiProjectors}
\alias{semiProjectors}
\title{Semi-Projectors}
\description{
 A function to compute the \eqn{L} and \eqn{R} semi-projectors for a given partition of the vertices.
}
\usage{
semiProjectors(gr, matrix = c("symmetric", "laplacian", "stochastic"), 
		norm = "row", sparse = TRUE, p = NULL)
}
\arguments{
  \item{gr}{A vector of \code{nrow(X)} or \code{vcount(X)} integers giving the group 
  label of every vertex in the partition.}

  \item{matrix}{The type of semi-projectors. For now \dQuote{symmetric}, \dQuote{laplacian} 
  and \dQuote{stochastic} are available.}
  
  \item{norm}{Either \dQuote{row} or \dQuote{col}. If set to \dQuote{row} the rows of the Laplacian matrix sum up to zero
  and the rows of the stochastic sum up to one; otherwise it is the columns.}
  
   \item{sparse}{When \code{X} is a graph, the \code{sparse} argument passed to \code{get.adjacency} of 
  the \pkg{igraph} package. Must be \code{TRUE} or \code{FALSE}.}

   \item{p}{A probability vector of length \code{length(gr)}. \code{p} is the stationary probability 
  distribution of a Markov chain when \code{matrix} = \dQuote{stochastic}.
  This parameter is ignored in all other cases.}
}

\details{
 The three types of semi-projectors are defined as follows.
 Let \eqn{\gamma(j)} label the group of vertex \eqn{j} in a partition of all the vertices.
 
The symmetric semi-projectors are defined as
 \deqn{L_{\alpha j}=R_{\alpha j}=\frac{1}{\sqrt{|\alpha|}}\delta_{\alpha\gamma(j)},}
 the (row) Laplacian semi-projectors as 
 \deqn{L_{\alpha j}=\frac{1}{|\alpha|}\delta_{\alpha\gamma(j)}\,\,\,\,
 \textmf{and}\,\,\,\, R_{\alpha j}=\delta_{\alpha\gamma(j)},}
 and the (row) stochastic semi-projectors as
 \deqn{L_{\alpha j}=\frac{p_{1}(j)}{\sum_{k\in\gamma(j)}p_{1}(k)}\delta_{\alpha\gamma(j)}\,\,\,\,
 \text{and}\,\,\,\, R_{\alpha j}=\delta_{\alpha\gamma(j)},}
where \eqn{p_1} is the (left) eigenvector associated with the one-eigenvalue of the stochastic matrix.
\eqn{L} and \eqn{R} are defined in a symmetric way when \code{norm = col}.
 All these semi-projectors verify various properties described in the reference. 
}

\value{
 \item{L}{The semi-projector \eqn{L}.}
 \item{R}{The semi-projector \eqn{R.}}
}

\references{
D. Morton de Lachapelle, D. Gfeller, and P. De Los Rios, Shrinking Matrices while Preserving their Eigenpairs with Application to the Spectral Coarse Graining of Graphs. Submitted to \emph{SIAM Journal on Matrix Analysis and Applications}, 2008.
\url{http://people.epfl.ch/david.morton}
}
\author{David Morton de Lachapelle \email{david.morton@epfl.ch}, \email{david.mortondelachapelle@swissquote.ch}}
\seealso{\code{\link{scg}}, \code{\link{normEps}}, \code{\link{grouping}}}

\examples{
# compute the semi-projectors and projector for the partition
# provided by a community detection method
g <- barabasi.game(20, m=1.5)
eb <- edge.betweenness.community(g)
memb <- community.to.membership(g, eb$merges, steps=15)$membership
lr <- semiProjectors(memb)
#In the symmetric case L = R
tcrossprod(lr$R)  #same as lr$R %*% t(lr$R)
P <- crossprod(lr$R)  #same as t(lr$R) %*% lr$R
#P is an orthogonal projector
isSymmetric(P)
sum( (P\%*\%P-P)^2 )

## use L and R to coarse-grain the graph Laplacian
lr <- semiProjectors(memb, matrix="laplacian")
L <- laplacianMatrix(g)
Lt <- lr$L\%*\%L\%*\%t(lr$R)  ##or better lr$L \%*\% tcrossprod(L,lr$R)
rowSums(Lt)
}
\keyword{array}
\keyword{graphs}

