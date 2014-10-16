#   IGraph R package
#   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
#   334 Harvard street, Cambridge, MA 02139 USA
#   
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
#   02110-1301 USA
#
###################################################################



#' Stochastic matrix of a graph
#' 
#' Retrieves the stochastic matrix of a graph of class \code{igraph}.
#' 
#' Let \eqn{M} be an \eqn{n \times n}{n x n} adjacency matrix with real
#' non-negative entries. Let us define \eqn{D = \textrm{diag}(\sum_{i}M_{1i},
#' \dots, \sum_{i}M_{ni})}{D=diag( sum(M[1,i], i), ..., sum(M[n,i], i) )}
#' 
#' The (row) stochastic matrix is defined as \deqn{W = D^{-1}M,}{W = inv(D) M,}
#' where it is assumed that \eqn{D} is non-singular.  Column stochastic
#' matrices are defined in a symmetric way.
#'
#' @aliases get.stochastic
#' @param graph The input graph. Must be of class \code{igraph}.
#' @param column.wise If \code{FALSE}, then the rows of the stochastic matrix
#' sum up to one; otherwise it is the columns.
#' @param sparse Logical scalar, whether to return a sparse matrix. The
#' \code{Matrix} package is needed for sparse matrices.
#' @return A regular matrix or a matrix of class \code{Matrix} if a
#' \code{sparse} argument was \code{TRUE}.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{as_adj}}
#' @export
#' @keywords graphs
#' @examples
#' 
#' library(Matrix)
#' ## g is a large sparse graph
#' g <- barabasi.game(n = 10^5, power = 2, directed = FALSE)
#' W <- stochastic_matrix(g, sparse=TRUE)
#' 
#' ## a dense matrix here would probably not fit in the memory
#' class(W)
#' 
#' ## may not be exactly 1, due to numerical errors
#' max(abs(rowSums(W))-1)
#' 
stochastic_matrix <- function(graph, column.wise=FALSE,
                           sparse=igraph_opt("sparsematrices")) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
 
  column.wise <- as.logical(column.wise)
  if (length(column.wise) != 1) {
    stop("`column.wise' must be a logical scalar")
  }

  sparse <- as.logical(sparse)
  if (length(sparse) != 1) {
    stop("`sparse' must be a logical scalar")
  }

  on.exit(.Call("R_igraph_finalizer", PACKAGE="igraph"))
  if (sparse) {
    res <- .Call("R_igraph_get_stochastic_sparsemat", graph, column.wise,
                 PACKAGE="igraph")
    res <- igraph.i.spMatrix(res)
  } else {
    res <- .Call("R_igraph_get_stochastic", graph, column.wise,
                 PACKAGE="igraph")
  }

  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    rownames(res) <- colnames(res) <- V(graph)$name
  }

  res
} 



#' SCG Problem Solver
#' 
#' This function solves the Spectral Coarse Graining (SCG) problem; either
#' exactly, or approximately but faster.
#' 
#' The algorithm \dQuote{optimum} solves exactly the SCG problem for each
#' eigenvector in \code{V}. The running time of this algorithm is \eqn{O(\max
#' nt \cdot m^2)}{O(max(nt) m^2)} for the symmetric and laplacian matrix
#' problems (i.e. when \code{mtype} is \dQuote{symmetric} or
#' \dQuote{laplacian}. It is \eqn{O(m^3)} for the stochastic problem. Here
#' \eqn{m} is the number of rows in \code{V}.  In all three cases, the memory
#' usage is \eqn{O(m^2)}.
#' 
#' The algorithms \dQuote{interv} and \dQuote{interv\_km} solve approximately
#' the SCG problem by performing a (for now) constant binning of the components
#' of the eigenvectors, that is \code{nt[i]} constant-size bins are used to
#' partition \code{V[,i]}. When \code{algo} = \dQuote{interv\_km}, the (Lloyd)
#' k-means algorithm is run on each partition obtained by \dQuote{interv} to
#' improve accuracy.
#' 
#' Once a minimizing partition (either exact or approximate) has been found for
#' each eigenvector, the final grouping is worked out as follows: two vertices
#' are grouped together in the final partition if they are grouped together in
#' each minimizing partition. In general the size of the final partition is not
#' known in advance when \code{ncol(V)}>1.
#' 
#' Finally, the algorithm \dQuote{exact\_scg} groups the vertices with equal
#' components in each eigenvector. The last three algorithms essentially have
#' linear running time and memory load.
#'
#' @aliases scgGrouping
#' @param V A numeric matrix of (eigen)vectors to be preserved by the coarse
#' graining (the vectors are to be stored column-wise in \code{V}).
#' @param nt A vector of positive integers of length one or equal to
#' \code{length(ev)}. When \code{algo} = \dQuote{optimum}, \code{nt} contains
#' the number of groups used to partition each eigenvector separately. When
#' \code{algo} is equal to \dQuote{interv\_km} or \dQuote{interv}, \code{nt}
#' contains the number of intervals used to partition each eigenvector. The
#' same partition size or number of intervals is used for each eigenvector if
#' \code{nt} is a single integer. When \code{algo} = \dQuote{exact\_cg} this
#' parameter is ignored.
#' @param mtype The type of semi-projectors used in the SCG. For now
#' \dQuote{symmetric}, \dQuote{laplacian} and \dQuote{stochastic} are
#' available.
#' @param algo The algorithm used to solve the SCG problem. Possible values are
#' \dQuote{optimum}, \dQuote{interv\_km}, \dQuote{interv} and
#' \dQuote{exact\_scg}.
#' @param p A probability vector of length equal to \code{nrow(V)}. \code{p} is
#' the stationary probability distribution of a Markov chain when \code{mtype}
#' = \dQuote{stochastic}. This parameter is ignored in all other cases.
#' @param maxiter A positive integer giving the maximum number of iterations of
#' the k-means algorithm when \code{algo} = \dQuote{interv\_km}. This parameter
#' is ignored in all other cases.
#' @return A vector of \code{nrow(V)} integers giving the group label of each
#' object (vertex) in the partition.
#' @author David Morton de Lachapelle \email{david.morton@@epfl.ch},
#' \email{david.mortondelachapelle@@swissquote.ch}
#' @seealso \link{SCG} for a detailed introduction. \code{\link{scg}},
#' \code{\link{scg_eps}}
#' @references D. Morton de Lachapelle, D. Gfeller, and P. De Los Rios,
#' Shrinking Matrices while Preserving their Eigenpairs with Application to the
#' Spectral Coarse Graining of Graphs. Submitted to \emph{SIAM Journal on
#' Matrix Analysis and Applications}, 2008.
#' \url{http://people.epfl.ch/david.morton}
#' @export
#' @keywords graphs
#' @examples
#' 
#' 
#' ## We are not running these examples any more, because they
#' ## take a long time to run and this is against the CRAN repository
#' ## policy. Copy and paste them by hand to your R prompt if
#' ## you want to run them.
#' 
#' \dontrun{
#' # eigenvectors of a random symmetric matrix
#' M <- matrix(rexp(10^6), 10^3, 10^3)
#' M <- (M + t(M))/2
#' V <- eigen(M, symmetric=TRUE)$vectors[,c(1,2)]
#' 
#' # displays size of the groups in the final partition
#' gr <- scg_group(V, nt=c(2,3))
#' col <- rainbow(max(gr))
#' plot(table(gr), col=col, main="Group size", xlab="group", ylab="size")
#' 
#' ## comparison with the grouping obtained by kmeans
#' ## for a partition of same size
#' gr.km <- kmeans(V,centers=max(gr), iter.max=100, nstart=100)$cluster
#' op <- par(mfrow=c(1,2))
#' plot(V[,1], V[,2], col=col[gr],
#' 	main = "SCG grouping",
#' 	xlab = "1st eigenvector",
#' 	ylab = "2nd eigenvector")
#' plot(V[,1], V[,2], col=col[gr.km],
#' 	main = "K-means grouping",
#' 	xlab = "1st eigenvector",
#' 	ylab = "2nd eigenvector")
#' par(op)
#' ## kmeans disregards the first eigenvector as it
#' ## spreads a much smaller range of values than the second one
#' 
#' ### comparing optimal and k-means solutions
#' ### in the one-dimensional case.
#' x <- rexp(2000, 2)
#' gr.true <- scg_group(cbind(x), 100)
#' gr.km <- kmeans(x, 100, 100, 300)$cluster
#' scg_eps(cbind(x), gr.true)
#' scg_eps(cbind(x), gr.km)
#' }
#' 
scg_group <- function(V, nt,
                         mtype=c("symmetric", "laplacian",
                           "stochastic"),
                         algo=c("optimum", "interv_km", "interv",
                           "exact_scg"),
                         p=NULL, maxiter=100) {

  V <- as.matrix(structure(as.double(V), dim=dim(V)))
  groups <- as.numeric(nt)

  mtype <- switch(igraph.match.arg(mtype), "symmetric"=1, 
                        "laplacian"=2, "stochastic"=3)
  algo <- switch(igraph.match.arg(algo), "optimum"=1,
                      "interv_km"=2, "interv"=3, "exact_scg"=4)
  if (!is.null(p)) p <- as.numeric(p)
  maxiter <- as.integer(maxiter)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_scg_grouping", V, as.integer(nt[1]),
               if (length(nt)==1) NULL else nt,
               mtype, algo, p, maxiter,
               PACKAGE="igraph")
  res
}



#' Semi-Projectors
#' 
#' A function to compute the \eqn{L} and \eqn{R} semi-projectors for a given
#' partition of the vertices.
#' 
#' The three types of semi-projectors are defined as follows.  Let
#' \eqn{\gamma(j)}{gamma(j)} label the group of vertex \eqn{j} in a partition
#' of all the vertices.
#' 
#' The symmetric semi-projectors are defined as \deqn{L_{\alpha j}=R_{\alpha
#' j}= }{% L[alpha,j] = R[alpha,j] = 1/sqrt(|alpha|)
#' delta[alpha,gamma(j)],}\deqn{
#' \frac{1}{\sqrt{|\alpha|}}\delta_{\alpha\gamma(j)},}{% L[alpha,j] =
#' R[alpha,j] = 1/sqrt(|alpha|) delta[alpha,gamma(j)],} the (row) Laplacian
#' semi-projectors as \deqn{L_{\alpha
#' j}=\frac{1}{|\alpha|}\delta_{\alpha\gamma(j)}\,\,\,\, }{% L[alpha,j] =
#' 1/|alpha| delta[alpha,gamma(j)] and R[alpha,j] =
#' delta[alpha,gamma(j)],}\deqn{ \textrm{and}\,\,\,\, R_{\alpha
#' j}=\delta_{\alpha\gamma(j)},}{% L[alpha,j] = 1/|alpha| delta[alpha,gamma(j)]
#' and R[alpha,j] = delta[alpha,gamma(j)],} and the (row) stochastic
#' semi-projectors as \deqn{L_{\alpha
#' j}=\frac{p_{1}(j)}{\sum_{k\in\gamma(j)}p_{1}(k)}\,\,\,\, }{% L[alpha,j] =
#' p[1][j] / sum(p[1][k]; k in gamma(j)) delta[alpha,gamma(j)] and R[alpha,j] =
#' delta[alpha,gamma(j)],}\deqn{ \textrm{and}\,\,\,\, R_{\alpha
#' j}=\delta_{\alpha\gamma(j)\delta_{\alpha\gamma(j)}},}{% L[alpha,j] = p[1][j]
#' / sum(p[1][k]; k in gamma(j)) delta[alpha,gamma(j)] and R[alpha,j] =
#' delta[alpha,gamma(j)],} where \eqn{p_1}{p[1]} is the (left) eigenvector
#' associated with the one-eigenvalue of the stochastic matrix. \eqn{L} and
#' \eqn{R} are defined in a symmetric way when \code{norm = col}. All these
#' semi-projectors verify various properties described in the reference.
#'
#' @aliases scgSemiProjectors
#' @param groups A vector of \code{nrow(X)} or \code{vcount(X)} integers giving
#' the group label of every vertex in the partition.
#' @param mtype The type of semi-projectors. For now \dQuote{symmetric},
#' \dQuote{laplacian} and \dQuote{stochastic} are available.
#' @param p A probability vector of length \code{length(gr)}. \code{p} is the
#' stationary probability distribution of a Markov chain when \code{mtype} =
#' \dQuote{stochastic}. This parameter is ignored in all other cases.
#' @param norm Either \dQuote{row} or \dQuote{col}. If set to \dQuote{row} the
#' rows of the Laplacian matrix sum up to zero and the rows of the stochastic
#' sum up to one; otherwise it is the columns.
#' @param sparse Logical scalar, whether to return sparse matrices.
#' @return \item{L}{The semi-projector \eqn{L}.} \item{R}{The semi-projector
#' \eqn{R}.}
#' @author David Morton de Lachapelle,
#' \url{http://people.epfl.ch/david.morton}.
#' @seealso \link{SCG} for a detailed introduction. \code{\link{scg}},
#' \code{\link{scg_eps}}, \code{\link{scg_group}}
#' @references D. Morton de Lachapelle, D. Gfeller, and P. De Los Rios,
#' Shrinking Matrices while Preserving their Eigenpairs with Application to the
#' Spectral Coarse Graining of Graphs. Submitted to \emph{SIAM Journal on
#' Matrix Analysis and Applications}, 2008.
#' \url{http://people.epfl.ch/david.morton}
#' @export
#' @examples
#' 
#' library(Matrix)
#' # compute the semi-projectors and projector for the partition
#' # provided by a community detection method
#' g <- barabasi.game(20, m=1.5)
#' eb <- cluster_edge_betweenness(g)
#' memb <- membership(eb)
#' lr <- scg_semi_proj(memb)
#' #In the symmetric case L = R
#' tcrossprod(lr$R)  # same as lr$R %*% t(lr$R)
#' P <- crossprod(lr$R)  # same as t(lr$R) %*% lr$R
#' #P is an orthogonal projector
#' isSymmetric(P)
#' sum( (P %*% P-P)^2 )
#' 
#' ## use L and R to coarse-grain the graph Laplacian
#' lr <- scg_semi_proj(memb, mtype="laplacian")
#' L <- laplacian_matrix(g)
#' Lt <- lr$L %*% L %*% t(lr$R)
#' ## or better lr$L %*% tcrossprod(L,lr$R)
#' rowSums(Lt)
#' 
scg_semi_proj <- function(groups,
                               mtype=c("symmetric", "laplacian",
                                 "stochastic"), p=NULL,
                               norm=c("row", "col"),
                               sparse=igraph_opt("sparsematrices")) {
  # Argument checks
  groups <- as.numeric(groups)-1
  mtype <- switch(igraph.match.arg(mtype), "symmetric"=1, 
  "laplacian"=2, "stochastic"=3)
  if (!is.null(p)) p <- as.numeric(p)
  norm <- switch(igraph.match.arg(norm), "row"=1, "col"=2)
  sparse <- as.logical(sparse)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_scg_semiprojectors", groups, mtype, p, norm,
               sparse,
               PACKAGE="igraph")

  if (sparse) {
    res$L <- igraph.i.spMatrix(res$L)
    res$R <- igraph.i.spMatrix(res$R)
  }
                
  res
}



#' All-in-one Function for the SCG of Matrices and Graphs
#' 
#' This function handles all the steps involved in the Spectral Coarse Graining
#' (SCG) of some matrices and graphs as described in the reference below.
#' 
#' Please see \link{SCG} for an introduction.
#' 
#' In the following \eqn{V} is the matrix of eigenvectors for which the SCG is
#' solved. \eqn{V} is calculated from \code{X}, if it is not given in the
#' \code{evec} argument.
#' 
#' The algorithm \dQuote{optimum} solves exactly the SCG problem for each
#' eigenvector in \code{V}. The running time of this algorithm is \eqn{O(\max
#' nt \cdot m^2)}{O(max(nt) m^2)} for the symmetric and laplacian matrix
#' problems (i.e. when \code{mtype} is \dQuote{symmetric} or
#' \dQuote{laplacian}. It is \eqn{O(m^3)} for the stochastic problem. Here
#' \eqn{m} is the number of rows in \code{V}.  In all three cases, the memory
#' usage is \eqn{O(m^2)}.
#' 
#' The algorithms \dQuote{interv} and \dQuote{interv\_km} solve approximately
#' the SCG problem by performing a (for now) constant binning of the components
#' of the eigenvectors, that is \code{nt[i]} constant-size bins are used to
#' partition \code{V[,i]}. When \code{algo} = \dQuote{interv\_km}, the (Lloyd)
#' k-means algorithm is run on each partition obtained by \dQuote{interv} to
#' improve accuracy.
#' 
#' Once a minimizing partition (either exact or approximate) has been found for
#' each eigenvector, the final grouping is worked out as follows: two vertices
#' are grouped together in the final partition if they are grouped together in
#' each minimizing partition. In general the size of the final partition is not
#' known in advance when \code{ncol(V)}>1.
#' 
#' Finally, the algorithm \dQuote{exact\_scg} groups the vertices with equal
#' components in each eigenvector. The last three algorithms essentially have
#' linear running time and memory load.
#' 
#' @param X The input graph or square matrix. Can be of class \code{igraph},
#' \code{matrix} or \code{Matrix}.
#' @param ev A vector of positive integers giving the indexes of the eigenpairs
#' to be preserved. For real eigenpairs, 1 designates the eigenvalue with
#' largest algebraic value, 2 the one with second largest algebraic value, etc.
#' In the complex case, it is the magnitude that matters.
#' @param nt A vector of positive integers of length one or equal to
#' \code{length(ev)}. When \code{algo} = \dQuote{optimum}, \code{nt} contains
#' the number of groups used to partition each eigenvector separately. When
#' \code{algo} is equal to \dQuote{interv\_km} or \dQuote{interv}, \code{nt}
#' contains the number of intervals used to partition each eigenvector. The
#' same partition size or number of intervals is used for each eigenvector if
#' \code{nt} is a single integer. When \code{algo} = \dQuote{exact\_cg} this
#' parameter is ignored.
#' @param groups A vector of \code{nrow(X)} or \code{vcount(X)} integers
#' labeling each group vertex in the partition. If this parameter is supplied
#' most part of the function is bypassed.
#' @param mtype Character scalar. The type of semi-projector to be used for the
#' SCG. For now \dQuote{symmetric}, \dQuote{laplacian} and \dQuote{stochastic}
#' are available.
#' @param algo Character scalar. The algorithm used to solve the SCG problem.
#' Possible values are \dQuote{optimum}, \dQuote{interv\_km}, \dQuote{interv}
#' and \dQuote{exact\_scg}.
#' @param norm Character scalar. Either \dQuote{row} or \dQuote{col}. If set to
#' \dQuote{row} the rows of the Laplacian matrix sum up to zero and the rows of
#' the stochastic matrix sum up to one; otherwise it is the columns.
#' @param direction Character scalar. When set to \dQuote{right}, resp.
#' \dQuote{left}, the parameters \code{ev} and \code{evec} refer to right,
#' resp. left eigenvectors. When passed \dQuote{default} it is the SCG
#' described in the reference below that is applied (common usage). This
#' argument is currently not implemented, and right eigenvectors are always
#' used.
#' @param evec A numeric matrix of (eigen)vectors to be preserved by the coarse
#' graining (the vectors are to be stored column-wise in \code{evec}). If
#' supplied, the eigenvectors should correspond to the indexes in \code{ev} as
#' no cross-check will be done.
#' @param p A probability vector of length \code{nrow(X)} (or
#' \code{vcount(X)}). \code{p} is the stationary probability distribution of a
#' Markov chain when \code{mtype} = \dQuote{stochastic}. This parameter is
#' ignored in all other cases.
#' @param use.arpack Logical scalar. When set to \code{TRUE} uses the function
#' \code{\link{arpack}} to compute eigenpairs. This parameter should be set to
#' \code{TRUE} if one deals with large (over a few thousands) AND sparse graphs
#' or matrices. This argument is not implemented currently and LAPACK is used
#' for solving the eigenproblems.
#' @param maxiter A positive integer giving the maximum number of iterations
#' for the k-means algorithm when \code{algo} = \dQuote{interv\_km}. This
#' parameter is ignored in all other cases.
#' @param sparse Logical scalar. Whether to return sparse matrices in the
#' result, if matrices are requested.
#' @param output Character scalar. Set this parameter to \dQuote{default} to
#' retrieve a coarse-grained object of the same class as \code{X}.
#' @param semproj Logical scalar. Set this parameter to \code{TRUE} to retrieve
#' the semi-projectors of the SCG.
#' @param epairs Logical scalar. Set this to \code{TRUE} to collect the
#' eigenpairs computed by \code{scg}.
#' @param stat.prob Logical scalar. This is to collect the stationary
#' probability \code{p} when dealing with stochastic matrices.
#' @return \item{Xt}{The coarse-grained graph, or matrix, possibly a sparse
#' matrix.} \item{groups}{A vector of \code{nrow(X)} or \code{vcount(X)}
#' integers giving the group label of each object (vertex) in the partition.}
#' \item{L}{The semi-projector \eqn{L} if \code{semproj = TRUE}.} \item{R}{The
#' semi-projector \eqn{R} if \code{semproj = TRUE}.} \item{values}{The computed
#' eigenvalues if \code{epairs = TRUE}.} \item{vectors}{The computed or
#' supplied eigenvectors if \code{epairs = TRUE}.} \item{p}{The stationary
#' probability vector if \code{mtype = stochastic} and \code{stat.prob = TRUE}.
#' For other matrix types this is missing.}
#' @author David Morton de Lachapelle,
#' \url{http://people.epfl.ch/david.morton}.
#' @seealso \link{SCG} for an introduction.  \code{\link{scg_eps}},
#' \code{\link{scg_group}} and \code{\link{scg_semi_proj}}.
#' @references D. Morton de Lachapelle, D. Gfeller, and P. De Los Rios,
#' Shrinking Matrices while Preserving their Eigenpairs with Application to the
#' Spectral Coarse Graining of Graphs. Submitted to \emph{SIAM Journal on
#' Matrix Analysis and Applications}, 2008.
#' \url{http://people.epfl.ch/david.morton}
#' @export
#' @keywords graphs
#' @examples
#' 
#' 
#' ## We are not running these examples any more, because they
#' ## take a long time (~20 seconds) to run and this is against the CRAN
#' ## repository policy. Copy and paste them by hand to your R prompt if
#' ## you want to run them.
#' 
#' \dontrun{
#' # SCG of a toy network
#' g <- make_full_graph(5) %du% make_full_graph(5) %du% make_full_graph(5)
#' g <- add_edges(g, c(1,6, 1,11, 6, 11))
#' cg <- scg(g, 1, 3, algo="exact_scg")
#' 
#' #plot the result
#' layout <- layout_with_kk(g)
#' nt <- vcount(cg$Xt)
#' col <- rainbow(nt)
#' vsize <- table(cg$groups)
#' ewidth <- round(E(cg$Xt)$weight,2)
#' 
#' op <- par(mfrow=c(1,2))
#' plot(g, vertex.color = col[cg$groups], vertex.size = 20,
#' 		vertex.label = NA, layout = layout)
#' plot(cg$Xt, edge.width = ewidth, edge.label = ewidth, 
#' 	vertex.color = col, vertex.size = 20*vsize/max(vsize),
#' 	vertex.label=NA, layout = layout_with_kk)
#' par(op)
#' 
#' ## SCG of real-world network
#' library(igraphdata)
#' data(immuno)
#' summary(immuno)
#' n <- vcount(immuno)
#' interv <- c(100,100,50,25,12,6,3,2,2)
#' cg <- scg(immuno, ev= n-(1:9), nt=interv, mtype="laplacian",
#'                         algo="interv", epairs=TRUE)
#' 
#' ## are the eigenvalues well-preserved?
#' gt <- cg$Xt
#' nt <- vcount(gt)
#' Lt <- laplacian_matrix(gt)
#' evalt <- eigen(Lt, only.values=TRUE)$values[nt-(1:9)]
#' res <- cbind(interv, cg$values, evalt)
#' res <- round(res,5)
#' colnames(res) <- c("interv","lambda_i","lambda_tilde_i")
#' rownames(res) <- c("N-1","N-2","N-3","N-4","N-5","N-6","N-7","N-8","N-9")
#' print(res)
#' 
#' ## use SCG to get the communities
#' com <- scg(laplacian_matrix(immuno), ev=n-c(1,2), nt=2)$groups
#' col <- rainbow(max(com))
#' layout <- layout_nicely(immuno)
#' 
#' plot(immuno, layout=layout, vertex.size=3, vertex.color=col[com],
#'                 vertex.label=NA)
#' 
#' ## display the coarse-grained graph
#' gt <- simplify(as.undirected(gt))
#' layout.cg <- layout_with_kk(gt)
#' com.cg <- scg(laplacian_matrix(gt), nt-c(1,2), 2)$groups
#' vsize <- sqrt(as.vector(table(cg$groups)))
#' 
#' op <- par(mfrow=c(1,2))
#' plot(immuno, layout=layout, vertex.size=3, vertex.color=col[com],
#'                 vertex.label=NA)
#' plot(gt, layout=layout.cg, vertex.size=15*vsize/max(vsize), 
#'                 vertex.color=col[com.cg],vertex.label=NA)
#' par(op)
#' 
#' }
#' 
scg <- function(X, ev, nt, groups=NULL, 
                mtype=c("symmetric", "laplacian", "stochastic"),
                algo=c("optimum", "interv_km", "interv",
                  "exact_scg"), norm=c("row", "col"),
                direction=c("default", "left", "right"),
                evec=NULL, p=NULL, use.arpack=FALSE, maxiter=300,
                sparse=igraph_opt("sparsematrices"),
                output=c("default", "matrix", "graph"), semproj=FALSE,
                epairs=FALSE, stat.prob=FALSE)
  UseMethod("scg")

#' @method scg igraph
#' @export

scg.igraph <- function(X, ev, nt, groups=NULL,
                       mtype=c("symmetric", "laplacian", "stochastic"),
                       algo=c("optimum", "interv_km", "interv",
                         "exact_scg"), norm=c("row", "col"),
                       direction=c("default", "left", "right"),
                       evec=NULL, p=NULL, use.arpack=FALSE, maxiter=300,
                       sparse=igraph_opt("sparsematrices"),
                       output=c("default", "matrix", "graph"), semproj=FALSE,
                       epairs=FALSE, stat.prob=FALSE) {
  
  myscg(graph=X, matrix=NULL, sparsemat=NULL, ev=ev, nt=nt,
        groups=groups, mtype=mtype, algo=algo,
        norm=norm, direction=direction, evec=evec, p=p,
        use.arpack=use.arpack, maxiter=maxiter, sparse=sparse,
        output=output, semproj=semproj, epairs=epairs,
        stat.prob=stat.prob)
}

#' @method scg matrix
#' @export

scg.matrix <- function(X, ev, nt, groups=NULL,
                       mtype=c("symmetric", "laplacian", "stochastic"),
                       algo=c("optimum", "interv_km", "interv",
                         "exact_scg"), norm=c("row", "col"),
                       direction=c("default", "left", "right"),
                       evec=NULL, p=NULL, use.arpack=FALSE, maxiter=300,
                       sparse=igraph_opt("sparsematrices"),
                       output=c("default", "matrix", "graph"), semproj=FALSE,
                       epairs=FALSE, stat.prob=FALSE) {
  
  myscg(graph=NULL, matrix=X, sparsemat=NULL, ev=ev, nt=nt,
        groups=groups, mtype=mtype, algo=algo,
        norm=norm, direction=direction, evec=evec, p=p, 
        use.arpack=use.arpack, maxiter=maxiter, sparse=sparse,
        output=output, semproj=semproj, epairs=epairs,
        stat.prob=stat.prob)
}

#' @method scg Matrix
#' @export

scg.Matrix <- function(X, ev, nt, groups=NULL,
                       mtype=c("symmetric", "laplacian", "stochastic"),
                       algo=c("optimum", "interv_km", "interv",
                         "exact_scg"), norm=c("row", "col"),
                       direction=c("default", "left", "right"),
                       evec=NULL, p=NULL, use.arpack=FALSE, maxiter=300,
                       sparse=igraph_opt("sparsematrices"),
                       output=c("default", "matrix", "graph"), semproj=FALSE,
                       epairs=FALSE, stat.prob=FALSE) {

  myscg(graph=NULL, matrix=NULL, sparsemat=X, ev=ev, nt=nt,
        groups=groups, mtype=mtype, algo=algo,
        norm=norm, direction=direction, evec=evec, p=p,
        use.arpack=use.arpack, maxiter=maxiter, sparse=sparse,
        output=output, semproj=semproj, epairs=epairs,
        stat.prob=stat.prob)
}

myscg <- function(graph, matrix, sparsemat, ev, nt, groups=NULL,
                  mtype=c("symmetric", "laplacian", "stochastic"),
                  algo=c("optimum", "interv_km", "interv",
                    "exact_scg"), norm=c("row", "col"),
                  direction=c("default", "left", "right"),
                  evec=NULL, p=NULL, use.arpack=FALSE, maxiter=300,
                  sparse=igraph_opt("sparsematrices"),
                  output=c("default", "matrix", "graph"), semproj=FALSE,
                  epairs=FALSE, stat.prob=FALSE) {

  ## Argument checks
  if (!is.null(graph))  { stopifnot(is_igraph(graph)) }
  if (!is.null(matrix)) { stopifnot(is.matrix(matrix)) }
  if (!is.null(sparsemat)) { stopifnot(inherits(sparsemat, "Matrix")) }

  if (!is.null(sparsemat)) { sparsemat <- as(sparsemat, "dgCMatrix") }
  ev <- as.numeric(as.integer(ev))
  nt <- as.numeric(as.integer(nt))
  if (!is.null(groups)) groups <- as.numeric(groups)
  mtype <- igraph.match.arg(mtype)
  algo <- switch(igraph.match.arg(algo), "optimum"=1,
                      "interv_km"=2, "interv"=3, "exact_scg"=4)
  if (!is.null(groups)) { storage.mode(groups) <- "double" }
  use.arpack <- as.logical(use.arpack)
  maxiter <- as.integer(maxiter)
  sparse <- as.logical(sparse)
  output <- switch(igraph.match.arg(output), "default"=1, "matrix"=2,
                   "graph"=3)
  semproj <- as.logical(semproj)
  epairs <- as.logical(epairs)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )

  if (mtype=="symmetric") {
    if (!is.null(evec)) { storage.mode(evec) <- "double" }
    res <- .Call("R_igraph_scg_adjacency", graph, matrix, sparsemat, ev,
                 nt, algo, evec, groups,
                 use.arpack, maxiter, sparse, output, semproj, epairs,
                 PACKAGE="igraph")

  } else if (mtype=="laplacian") {
    norm <- switch(igraph.match.arg(norm), "row"=1, "col"=2)
    if (!is.null(evec)) { storage.mode(evec) <- "complex" }
    direction <- switch(igraph.match.arg(direction), "default"=1, "left"=2,
                        "right"=3)
    res <- .Call("R_igraph_scg_laplacian", graph, matrix, sparsemat, ev,
                 nt, algo, norm, direction,
                 evec, groups, use.arpack, maxiter, sparse, output,
                 semproj, epairs,
                 PACKAGE="igraph")

  } else if (mtype=="stochastic") {
    norm <- switch(igraph.match.arg(norm), "row"=1, "col"=2)
    if (!is.null(evec)) { storage.mode(evec) <- "complex" }
    if (!is.null(p)) { storage.mode(p) <- "double" }
    stat.prob <- as.logical(stat.prob)
    res <- .Call("R_igraph_scg_stochastic", graph, matrix, sparsemat, ev,
                 nt, algo, norm, evec, groups, p, use.arpack,
                 maxiter, sparse, output, semproj, epairs, stat.prob,
                 PACKAGE="igraph")    
  }

  if (!is.null(res$Xt) &&
      class(res$Xt) == "igraph.tmp.sparse") {
    res$Xt <- igraph.i.spMatrix(res$Xt)
  }
  if (!is.null(res$L) && class(res$L) == "igraph.tmp.sparse") {
    res$L <- igraph.i.spMatrix(res$L)
  }
  if (!is.null(res$R) && class(res$R) == "igraph.tmp.sparse") {
    res$R <- igraph.i.spMatrix(res$R)
  }

  res
}
