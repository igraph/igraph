\name{scgGrouping}
\alias{scgGrouping}
\title{SCG Problem Solver}
\description{
  This function solves the Spectral Coarse Graining (SCG) problem;
  either exactly, or approximately but faster.
}
\usage{
scgGrouping(V, nt, mtype = c("symmetric", "laplacian",
          "stochastic"), algo = c("optimum", "interv_km",
          "interv","exact_scg"), p = NULL, maxiter = 100)
}
\arguments{
  \item{V}{A numeric matrix of (eigen)vectors to be preserved by the
    coarse graining (the vectors are to be stored column-wise in
    \code{V}).} 
  \item{nt}{A vector of positive integers of length one or equal to
    \code{length(ev)}. When \code{algo} = \dQuote{optimum}, \code{nt}
    contains the number of groups used to partition each eigenvector
    separately. When \code{algo} is equal to \dQuote{interv\_km} or
    \dQuote{interv}, \code{nt} contains the number of intervals used to
    partition each eigenvector. The same partition size or number of
    intervals is used for each eigenvector if \code{nt} is a single
    integer. When \code{algo} = \dQuote{exact\_cg} this parameter is
    ignored, but a value between 2 and vcount(X) (or nrow(X)) must be
    provided TODO: improve this.}
  \item{mtype}{The type of semi-projectors used in the SCG. For now
    \dQuote{symmetric}, \dQuote{laplacian} and \dQuote{stochastic} are
    available.}
  \item{algo}{The algorithm used to solve the SCG problem. Possible
    values are \dQuote{optimum}, \dQuote{interv\_km}, \dQuote{interv}
    and \dQuote{exact\_scg}.} 
  \item{p}{A probability vector of length equal to
    \code{nrow(V)}. \code{p} is the stationary probability distribution
    of a Markov chain when \code{mtype} = \dQuote{stochastic}. This
    parameter is ignored in all other cases.} \item{maxiter}{A positive
    integer giving the maximum number of iterations of the k-means
    algorithm when \code{algo} = \dQuote{interv\_km}. This parameter is
    ignored in all other cases.} 
}
\details{
  The algorithm \dQuote{optimum} solves exactly the SCG problem for each
  eigenvector in \code{V}. The running time of this algorithm is
  \eqn{O(\max nt \cdot m^2)}{O(max(nt) m^2)} for the symmetric and
  laplacian matrix problems (i.e. when \code{mtype} is
  \dQuote{symmetric} or \dQuote{laplacian}. It is \eqn{O(m^3)} for the
  stochastic problem. Here \eqn{m} is the number of rows in \code{V}.
  In all three cases, the memory usage is \eqn{O(m^2)}.

  The algorithms \dQuote{interv} and \dQuote{interv\_km} solve
  approximately the SCG problem by performing a (for now) constant
  binning of the components of the eigenvectors, that is \code{nt[i]}
  constant-size bins are used to partition \code{V[,i]}. When
  \code{algo} = \dQuote{interv\_km}, the (Lloyd) k-means algorithm is
  run on each partition obtained by \dQuote{interv} to improve
  accuracy. 

  Once a minimizing partition (either exact or approximate) has been
  found for each eigenvector, the final grouping is worked out as
  follows: two vertices are grouped together in the final partition if
  they are grouped together in each minimizing partition. In general the
  size of the final partition is not known in advance when
  \code{ncol(V)}>1. 

  Finally, the algorithm \dQuote{exact\_scg} groups the vertices with
  equal components in each eigenvector. The last three algorithms
  essentially have linear running time and memory load.
}
\value{
  A vector of \code{nrow(V)} integers giving the group label of each
  object (vertex) in the partition.
}
\references{
D. Morton de Lachapelle, D. Gfeller, and P. De Los Rios, Shrinking
  Matrices while Preserving their Eigenpairs with Application to the
  Spectral Coarse Graining of Graphs. Submitted to \emph{SIAM Journal on
    Matrix Analysis and Applications}, 2008.
  \url{http://people.epfl.ch/david.morton}
}
\author{David Morton de Lachapelle \email{david.morton@epfl.ch},
  \email{david.mortondelachapelle@swissquote.ch}} 
\seealso{\link{SCG} for a detailed introduction. \code{\link{scg}},
  \code{\link{scgNormEps}}}
\examples{
# eigenvectors of a random symmetric matrix
M <- matrix(rexp(10^6), 10^3, 10^3)
M <- (M + t(M))/2
V <- eigen(M, symmetric=TRUE)$vectors[,c(1,2)]

# displays size of the groups in the final partition
gr <- scgGrouping(V, nt=c(2,3))
col <- rainbow(max(gr))
plot(table(gr), col=col, main="Group size", xlab="group", ylab="size")

## comparison with the grouping obtained by kmeans
## for a partition of same size
gr.km <- kmeans(V,centers=max(gr), iter.max=100, nstart=100)$cluster
op <- par(mfrow=c(1,2))
plot(V[,1], V[,2], col=col[gr],
	main = "SCG grouping",
	xlab = "1st eigenvector",
	ylab = "2nd eigenvector")
plot(V[,1], V[,2], col=col[gr.km],
	main = "K-means grouping",
	xlab = "1st eigenvector",
	ylab = "2nd eigenvector")
par(op)
## kmeans disregards the first eigenvector as it
## spreads a much smaller range of values than the second one

### comparing optimal and k-means solutions
### in the one-dimensional case.
x <- rexp(2000, 2)
gr.true <- scgGrouping(cbind(x), 100)
gr.km <- kmeans(x, 100, 100, 300)$cluster
scgNormEps(cbind(x), gr.true)
scgNormEps(cbind(x), gr.km)
}
\keyword{graphs}
