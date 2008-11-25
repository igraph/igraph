\name{scg}
\alias{scg}
\title{All-in-one Function for the SCG of Matrices and Graphs}
\description{
This function handles all the steps involved in the Spectral Coarse Graining (SCG) of some matrices 
and graphs as described in the reference below.
}
\usage{
scg(X, ev, nt, matrix = c("symmetric", "laplacian", "stochastic"),
    algo = c("optimum", "interv_km", "interv", "exact_scg"),
    norm = "row", sparse = TRUE,
    direction = c("default", "left", "right"), evec = NULL,
    p = NULL, gr = NULL, use.arpack = FALSE, maxiter = 100, 
    output = c("default", "matrix", "graph"), semproj = FALSE, 
    epairs = FALSE, stat.prob = FALSE, ...)
}
\arguments{
  \item{X}{The input graph or square matrix. Can be of class \code{igraph}, \code{matrix} or \code{Matrix}.}
  
  \item{ev}{A vector of positive integers giving the indexes of the eigenpairs to be preserved.
  For real eigenpairs, 1 designates the eigenvalue with largest algebraic value, 2 the one with second largest algebraic value, etc.
  In the complex case, it is the magnitude that matters.}
  
  \item{nt}{A vector of positive integers of length one or equal to \code{length(ev)}.
  When \code{algo} = \dQuote{optimum}, \code{nt} contains the number of groups used to partition each eigenvector separately.
  When \code{algo} is equal to \dQuote{interv\_km} or \dQuote{interv}, \code{nt} contains the number of intervals
  used to partition each eigenvector.
  The same partition size or number of intervals is used for each eigenvector if \code{nt} is a single integer.
  When \code{algo} = \dQuote{exact\_cg} this parameter is ignored, but a value between 2 and \vcount(X) (or \nrow(X)) must be provided TODO: improve this.}
  
  \item{matrix}{The type of semi-projector to be used for the SCG. For now \dQuote{symmetric}, \dQuote{laplacian} 
  and \dQuote{stochastic} are available.}
  
  \item{algo}{The algorithm used to solve the SCG problem. Possible values are \dQuote{optimum}, \dQuote{interv\_km},
   \dQuote{interv} and \dQuote{exact\_scg}.}
  
  \item{norm}{Either \dQuote{row} or \dQuote{col}. If set to \dQuote{row} the rows of the Laplacian matrix sum up to zero
  and the rows of the stochastic matrix sum up to one; otherwise it is the columns.}
  
  \item{sparse}{When \code{X} is a graph, the \code{sparse} argument passed to \code{get.adjacency} of 
  the \pkg{igraph} package. Must be \code{TRUE} or \code{FALSE}.}
  
  \item{direction}{When set to \dQuote{right}, resp. \dQuote{left}, the parameters \code{ev} and \code{evec} refer to right, 
  resp. left eigenvectors. When passed \dQuote{default} it is the SCG described in the reference below that is applied (common usage).}
  
  \item{evec}{A numeric matrix of (eigen)vectors to be preserved by the coarse graining 
  (the vectors are to be stored column-wise in \code{evec}). If supplied, the eigenvectors should correspond to the indexes
   in \code{ev} as no cross-check will be done.}
  
  \item{p}{A probability vector of length \code{nrow(V)}. \code{p} is the stationary probability 
  distribution of a Markov chain when \code{matrix} = \dQuote{stochastic}.
  This parameter is ignored in all other cases.}
  
 \item{gr}{A vector of \code{nrow(X)} or \code{vcount(X)} integers labeling each group vertex in the partition. If this parameter is supplied most part of the function is bypassed.} 
  
  \item{use.arpack}{When set to \code{TRUE} uses the function \code{arpack} of the \pkg{igraph} package to compute eigenpairs. 
  This parameter should be set to \code{TRUE} if one deals with large (over a few thousands) AND sparse graphs or matrices.}
  
  \item{maxiter}{A positive integer giving the maximum number of iterations for the k-means algorithm 
  when \code{algo} = \dQuote{interv\_km}. This parameter is ignored in all other cases.}
  
  \item{output}{Set this parameter to \dQuote{default} to retrieve a coarse-grained object of the same class as \code{X}.}
  
  \item{semproj}{Set this parameter to \code{TRUE} to retrieve the semi-projectors of the SCG.}
  
  \item{epairs}{Set this to \code{TRUE} to collect the eigenpairs computed by \code{scg}.}
  
  \item{stat.prob}{This is to collect the stationary probability \code{p} when dealing with stochastic matrices.}
  
  \item{\dots}{Other arguments passed on to \code{get.adjacency} (only when \code{X} is a graph).}
}

\details{
The algorithm \dQuote{optimum} solves exactly the SCG problem for each eigenvector in \code{V}. The running time of this algorithm is \eqn{O(\code{nt[i]} \cdot \code{nrow(V)}^2)} when \code{matrix} is \dQuote{symmetric} or \dQuote{laplacian} and \eqn{O(\code{nrow(V)}^3)} when \code{matrix} = \dQuote{stochastic}. In all these cases, the memory usage is \eqn{O(\code{nrow(V)}^2)}.

The algorithms \dQuote{interv} and \dQuote{interv\_km} solve approximately the SCG problem by performing a (for now) constant binning of the components of the eigenvectors, that is \code{nt[i]} constant-size bins are used to partition \code{V[,i]}. When \code{algo} = \dQuote{interv\_km}, the (Lloyd) k-means algorithm is run on each partition obtained by \dQuote{interv} to improve accuracy.

Once a minimizing partition (either exact or approximate) has been found for each eigenvector, the final grouping is worked out as follows: two vertices are grouped together in the final partition if they are grouped together in each minimizing partition. In general the size of the final partition is not known in advance when \code{ncol(V)}>1.

Finally, the algorithm \dQuote{exact\_scg} groups the vertices with equal components in each eigenvector.
The last three algorithms essentially have linear running time and memory load.
}

\value{
  \item{Xt }{The coarse-grained graph or matrix, depending on  the class of \code{X} and the value of \code{output}.}
  \item{gr}{A vector of \code{nrow(X)} or \code{vcount(X)} integers giving the group label of each object (vertex)
 in the partition.}
 \item{L}{The semi-projector \eqn{L} if \code{semproj = TRUE}.}
 \item{R}{The semi-projector \eqn{R} if \code{semproj = TRUE}.}
 \item{values}{The computed eigenvalues if \code{epairs = TRUE}.}
 \item{vectors}{The computed or supplied eigenvectors if \code{epairs = TRUE}.}
 \item{p}{The stationary probability vector if \code{matrix = stochastic} and \code{stat.prob = TRUE}.}
}

\references{
D. Morton de Lachapelle, D. Gfeller, and P. De Los Rios, Shrinking Matrices while Preserving their Eigenpairs with Application to the Spectral Coarse Graining of Graphs. Submitted to \emph{SIAM Journal on Matrix Analysis and Applications}, 2008.
\url{http://people.epfl.ch/david.morton}
}
\author{David Morton de Lachapelle \email{david.morton@epfl.ch}, \email{david.mortondelachapelle@swissquote.ch}}
\seealso{\code{\link{scg}}, \code{\link{normEps}}, \code{\link{groupComposition}}, \code{\link{renumber}} }

\examples{
# SCG of a toy network
g <- graph.full(5) %du% graph.full(5) %du% graph.full(5)
g <- add.edges(g, c(0,5, 0,10, 5, 10))
cg <- scg(g, 1, 3)

#plot the result
layout <- layout.kamada.kawai(g)
nt <- vcount(cg$Xt)
col <- rainbow(nt)
vsize <- table(cg$gr)
ewidth <- round(E(cg$Xt)$weight,2)
op <- par(mfrow=c(1,2))
plot(g, vertex.color = col[cg$g], vertex.size = 20,
		vertex.label = NA, layout = layout)
plot(cg$Xt, edge.width = ewidth, edge.label = ewidth, 
	vertex.color = col, vertex.size = 20*vsize/max(vsize),
	vertex.label=NA, layout = layout.kamada.kawai)
par(op)

## SCG of real-world network
data(immuno)
summary(immuno)
n <- vcount(immuno)
interv <- c(100,100,50,25,12,6,3,2,2)
cg <- scg(immuno, ev= n-(1:9), nt=interv, matrix="laplacian",
			algo="interv", epairs=TRUE)
##TODO: clear Matrix warnings

## are the eigenvalues well-preserved?
gt <- cg$Xt
nt <- vcount(gt)
Lt <- laplacianMatrix(gt, attr="weight")
evalt <- eigen(Lt, only.values=TRUE)$values[nt-(1:9)]
res <- cbind(interv, cg$values, evalt)
res <- round(res,5)
colnames(res) <- c("interv","lambda_i","lambda_tilde_i")
rownames(res) <- c("N-1","N-2","N-3","N-4","N-5","N-6","N-7","N-8","N-9")
print(res)

## use SCG to get the communities
com <- scg(laplacianMatrix(immuno), n-c(1,2), 2)$gr
col <- rainbow(max(com))
layout <- layout.kamada.kawai(immuno) ## this can take some time
plot(immuno, layout=layout, vertex.size=3, vertex.color=col[com],
		vertex.label=NA)

## display the coarse-grained graph
gt <- simplify(as.undirected(gt))
layout.cg <- layout.kamada.kawai(gt)
com.cg <- scg(laplacianMatrix(gt), nt-c(1,2), 2)$gr
vsize <- sqrt(as.vector(table(cg$gr)))
op <- par(mfrow=c(1,2))
plot(immuno, layout=layout, vertex.size=3, vertex.color=col[com],
		vertex.label=NA)
plot(gt, layout=layout.cg, vertex.size=15*vsize/max(vsize), 
		vertex.color=col[com.cg],vertex.label=NA)
par(op)
}
\keyword{array}
\keyword{graphs}
