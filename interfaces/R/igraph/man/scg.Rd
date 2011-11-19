\name{scg}
\alias{scg}
\concept{Spectral coarse graining}
\title{All-in-one Function for the SCG of Matrices and Graphs}
\description{
  This function handles all the steps involved in the Spectral Coarse
  Graining (SCG) of some matrices and graphs as described in the
  reference below.
}
\usage{
scg(X, ev, nt, groups = NULL, mtype = c("symmetric", "laplacian",
      "stochastic"), algo = c("optimum", "interv_km", "interv",
      "exact_scg"), norm = c("row", "col"), direction = c("default",
      "left", "right"), evec = NULL, p = NULL, use.arpack = FALSE,
    maxiter = 300, sparse = getIgraphOpt("sparsematrices"), output =
    c("default", "matrix", "graph"), semproj = FALSE, epairs = FALSE,
    stat.prob = FALSE) 
}
\arguments{
  \item{X}{The input graph or square matrix. Can be of class
    \code{igraph}, \code{matrix} or \code{Matrix}.} 
  \item{ev}{A vector of positive integers giving the indexes of the
    eigenpairs to be preserved. For real eigenpairs, 1 designates the
    eigenvalue with largest algebraic value, 2 the one with second
    largest algebraic value, etc. In the complex case, it is the
    magnitude that matters.}
  \item{nt}{A vector of positive integers of length one or equal
    to \code{length(ev)}. When \code{algo} = \dQuote{optimum},
    \code{nt} contains the number of groups used to partition
    each eigenvector separately. When \code{algo} is equal to
    \dQuote{interv\_km} or \dQuote{interv}, \code{nt} contains
    the number of intervals used to partition each eigenvector. The same
    partition size or number of intervals is used for each eigenvector
    if \code{nt} is a single integer. When \code{algo} =
    \dQuote{exact\_cg} this parameter is ignored.}
  \item{groups}{A vector of \code{nrow(X)} or \code{vcount(X)} integers
    labeling each group vertex in the partition. If this parameter is
    supplied most part of the function is bypassed.}
  \item{mtype}{Character scalar. The type of semi-projector to be
    used for the SCG. For now \dQuote{symmetric}, \dQuote{laplacian} and
    \dQuote{stochastic} are available.}
  \item{algo}{Character scalar. The algorithm used to solve the SCG
    problem. Possible values are \dQuote{optimum}, \dQuote{interv\_km},
    \dQuote{interv} and \dQuote{exact\_scg}.}
  \item{norm}{Character scalar. Either \dQuote{row} or \dQuote{col}. If
    set to \dQuote{row} the rows of the Laplacian matrix sum up to zero
    and the rows of the stochastic matrix sum up to one; otherwise it is
    the columns.}  
  \item{direction}{Character scalar. When set to \dQuote{right},
    resp. \dQuote{left}, the parameters \code{ev} and \code{evec} refer
    to right, resp. left eigenvectors. When passed \dQuote{default} it
    is the SCG described in the reference below that is applied (common
    usage). This argument is currently not implemented, and right
    eigenvectors are always used.}
  \item{evec}{A numeric matrix of (eigen)vectors to be preserved by the
    coarse graining (the vectors are to be stored column-wise in
    \code{evec}). If supplied, the eigenvectors should correspond to the
    indexes in \code{ev} as no cross-check will be done.}
  \item{p}{A probability vector of length \code{nrow(X)} (or
    \code{vcount(X)}). \code{p} is the stationary probability
    distribution of a Markov chain when \code{mtype} =
    \dQuote{stochastic}. This parameter is ignored in all other cases.}
 \item{use.arpack}{Logical scalar. When set to \code{TRUE} uses the
   function \code{arpack} of the \pkg{igraph} package to compute
   eigenpairs. This parameter should be set to \code{TRUE} if one deals
   with large (over a few thousands) AND sparse graphs or matrices. This
   argument is not implemented currently and LAPACK is used for solving
   the eigenproblems. }
  \item{maxiter}{A positive integer giving the maximum number of
    iterations for the k-means algorithm when \code{algo} =
    \dQuote{interv\_km}. This parameter is ignored in all other cases.}
  \item{sparse}{Logical scalar. Whether to return sparse matrices in the
    result, if matrices are requested.}
  \item{output}{Character scalar. Set this parameter to \dQuote{default}
    to retrieve a coarse-grained object of the same class as \code{X}.}
  \item{semproj}{Logical scalar. Set this parameter to \code{TRUE} to
    retrieve the semi-projectors of the SCG.}
  \item{epairs}{Logical scalar. Set this to \code{TRUE} to collect the
    eigenpairs computed by \code{scg}.}
  \item{stat.prob}{Logical scalar. This is to collect the stationary
    probability \code{p} when dealing with stochastic matrices.}
}
\details{
  Please see \link{SCG} for an introduction.
  
  In the following \eqn{V} is the matrix of eigenvectors for which the
  SCG is solved. \eqn{V} is calculated from \code{X}, if it is not given
  in the \code{evec} argument.
  
  The algorithm \dQuote{optimum} solves exactly the SCG problem for each
  eigenvector in \code{V}. The running time of this algorithm is
  \eqn{O(\code{nt[i]} \cdot \code{nrow(V)}^2)} when
  \code{mtype} is \dQuote{symmetric} or \dQuote{laplacian} and
  \eqn{O(\code{nrow(V)}^3)} when \code{mtype} = \dQuote{stochastic}. In
  all these cases, the memory usage is \eqn{O(\code{nrow(V)}^2)}.

  The algorithms \dQuote{interv} and \dQuote{interv\_km} solve
  approximately the SCG problem by performing a (for now) constant
  binning of the components of the eigenvectors, that is
  \code{nt[i]} constant-size bins are used to partition
  \code{V[,i]}. When \code{algo} = \dQuote{interv\_km}, the (Lloyd)
  k-means algorithm is run on each partition obtained by \dQuote{interv}
  to improve accuracy. 

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
  \item{Xt}{The coarse-grained graph, or matrix, possibly a sparse matrix.}
  \item{groups}{A vector of \code{nrow(X)} or \code{vcount(X)} integers
    giving the group label of each object (vertex) in the partition.}
  \item{L}{The semi-projector \eqn{L} if \code{semproj = TRUE}.}
  \item{R}{The semi-projector \eqn{R} if \code{semproj = TRUE}.}
  \item{values}{The computed eigenvalues if \code{epairs = TRUE}.}
  \item{vectors}{The computed or supplied eigenvectors if \code{epairs =
      TRUE}.}
  \item{p}{The stationary probability vector if \code{mtype =
      stochastic} and \code{stat.prob = TRUE}. For other matrix types
    this is missing.}
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
\seealso{\link{SCG} for an introduction.
  \code{\link{scgNormEps}}, \code{\link{scgGrouping}} and
  \code{\link{scgSemiProjectors}}.}
\examples{
# SCG of a toy network
g <- graph.full(5) \%du\% graph.full(5) \%du\% graph.full(5)
g <- add.edges(g, c(1,6, 1,11, 6, 11))
cg <- scg(g, 1, 3)

#plot the result
layout <- layout.kamada.kawai(g)
nt <- vcount(cg$scg_graph)
col <- rainbow(nt)
vsize <- table(cg$groups)
ewidth <- round(E(cg$scg_graph)$weight,2)
\dontrun{
op <- par(mfrow=c(1,2))
plot(g, vertex.color = col[cg$groups], vertex.size = 20,
		vertex.label = NA, layout = layout)
plot(cg$scg_graph, edge.width = ewidth, edge.label = ewidth, 
	vertex.color = col, vertex.size = 20*vsize/max(vsize),
	vertex.label=NA, layout = layout.kamada.kawai)
par(op)
}

## SCG of real-world network
## TODO: add immuno example
}
\keyword{graphs}
