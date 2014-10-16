
## ----------------------------------------------------------------
##
##   IGraph R package
##   Copyright (C) 2005-2014  Gabor Csardi <csardi.gabor@gmail.com>
##   334 Harvard street, Cambridge, MA 02139 USA
##
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
##   02110-1301 USA
##
## -----------------------------------------------------------------

graph.adjacency.dense <- function(adjmatrix, mode=c("directed", "undirected", "max",
                                               "min", "upper", "lower", "plus"),
                                  weighted=NULL, diag=TRUE) {

  mode <- igraph.match.arg(mode)
  mode <- switch(mode,
                 "directed"=0,
                 "undirected"=1,
                 "max"=1,
                 "upper"=2,
                 "lower"=3,
                 "min"=4,
                 "plus"=5)

  mode(adjmatrix) <- "double"

  if (!is.null(weighted)) {
    if (is.logical(weighted) && weighted) {
      weighted <- "weight"
    }
    if (!is.character(weighted)) {
      stop("invalid value supplied for `weighted' argument, please see docs.")
    }

    if (nrow(adjmatrix) != ncol(adjmatrix)) {
      stop("not a square matrix")
    }

    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    res <- .Call("R_igraph_weighted_adjacency", adjmatrix,
                 as.numeric(mode), weighted, diag,
                 PACKAGE="igraph")
  } else {

    adjmatrix <- as.matrix(adjmatrix)
    attrs <- attributes(adjmatrix)
    adjmatrix <- as.numeric(adjmatrix)
    attributes(adjmatrix) <- attrs

    if (!diag) { diag(adjmatrix) <- 0 }

    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    res <- .Call("R_igraph_graph_adjacency", adjmatrix, as.numeric(mode),
                 PACKAGE="igraph")
  }

  res
}

graph.adjacency.sparse <- function(adjmatrix, mode=c("directed", "undirected", "max",
                                                "min", "upper", "lower", "plus"),
                                   weighted=NULL, diag=TRUE) {

  mode <- igraph.match.arg(mode)

  if (!is.null(weighted)) {
    if (is.logical(weighted) && weighted) {
      weighted <- "weight"
    }
    if (!is.character(weighted)) {
      stop("invalid value supplied for `weighted' argument, please see docs.")
    }
  }

  mysummary <- Matrix::summary

  if (nrow(adjmatrix) != ncol(adjmatrix)) {
    stop("not a square matrix")
  }

  vc <- nrow(adjmatrix)

  ## to remove non-redundancies that can persist in a dgtMatrix
  if(inherits(adjmatrix, "dgTMatrix")) {
    adjmatrix = as(adjmatrix, "CsparseMatrix")
  }

  if (is.null(weighted) && mode=="undirected") { mode <- "max" }

  if (mode == "directed") {
    ## DIRECTED
    el <- mysummary(adjmatrix)
    if (!diag) { el <- el[ el[,1] != el[,2], ] }
  } else if (mode == "undirected") {
    ## UNDIRECTED, must be symmetric if weighted
    if (!is.null(weighted) && !Matrix::isSymmetric(adjmatrix)) {
      stop("Please supply a symmetric matrix if you want to create a weighted graph with mode=UNDIRECTED.")
    }
    if (diag) {
      adjmatrix <- Matrix::tril(adjmatrix)
    } else {
      adjmatrix <- Matrix::tril(adjmatrix, -1)
    }
    el <- mysummary(adjmatrix)
  } else if (mode=="max") {
    ## MAXIMUM
    el <- mysummary(adjmatrix)
    rm(adjmatrix)
    if (!diag) { el <- el[ el[,1] != el[,2], ] }
    el <- el[ el[,3] != 0, ]
    w <- el[,3]
    el <- el[,1:2]
    el <- cbind( pmin(el[,1],el[,2]), pmax(el[,1], el[,2]) )
    o <- order(el[,1], el[,2])
    el <- el[o,,drop=FALSE]
    w <- w[o]
    if (nrow(el) > 1) {
      dd <- el[2:nrow(el),1] == el[1:(nrow(el)-1),1] &
        el[2:nrow(el),2] == el[1:(nrow(el)-1),2]
      dd <- which(dd)
      if (length(dd)>0) {
        mw <- pmax(w[dd], w[dd+1])
        w[dd] <- mw
        w[dd+1] <- mw
        el <- el[-dd,,drop=FALSE]
        w <- w[-dd]
      }
    }
    el <- cbind(el, w)
  } else if (mode=="upper") {
    ## UPPER
    if (diag) {
      adjmatrix <- Matrix::triu(adjmatrix)
    } else {
      adjmatrix <- Matrix::triu(adjmatrix, 1)
    }
    el <- mysummary(adjmatrix)
    rm(adjmatrix)
    if (!diag) { el <- el[ el[,1] != el[,2], ] }
  } else if (mode=="lower") {
    ## LOWER
    if (diag) {
      adjmatrix <- Matrix::tril(adjmatrix)
    } else {
      adjmatrix <- Matrix::tril(adjmatrix, -1)
    }
    el <- mysummary(adjmatrix)
    rm(adjmatrix)
    if (!diag) { el <- el[ el[,1] != el[,2], ] }
  } else if (mode=="min") {
    ## MINIMUM
    adjmatrix <- sign(adjmatrix) * sign(Matrix::t(adjmatrix)) * adjmatrix
    el <- mysummary(adjmatrix)
    if (!diag) { el <- el[ el[,1] != el[,2], ] }
    el <- el[ el[,3] != 0, ]
    w <- el[,3]
    el <- el[,1:2]
    el <- cbind( pmin(el[,1],el[,2]), pmax(el[,1], el[,2]) )
    o <- order(el[,1], el[,2])
    el <- el[o,]
    w <- w[o]
    if (nrow(el) > 1) {
      dd <- el[2:nrow(el),1] == el[1:(nrow(el)-1),1] &
        el[2:nrow(el),2] == el[1:(nrow(el)-1),2]
      dd <- which(dd)
      if (length(dd)>0) {
        mw <- pmin(w[dd], w[dd+1])
        w[dd] <- mw
        w[dd+1] <- mw
        el <- el[-dd,]
        w <- w[-dd]
      }
    }
    el <- cbind(el, w)
  } else if (mode=="plus") {
    ## PLUS
    adjmatrix <- adjmatrix + Matrix::t(adjmatrix)
    if (diag) {
      adjmatrix <- Matrix::tril(adjmatrix)
    } else {
      adjmatrix <- Matrix::tril(adjmatrix, -1)
    }
    el <- mysummary(adjmatrix)
    if (diag) {
      loop <- el[,1] == el[,2]
      el[loop,3] <- el[loop,3] / 2
    }
    el <- el[ el[,3] != 0, ]
    rm(adjmatrix)
  }

  if (!is.null(weighted)) {
    res <- make_empty_graph(n=vc, directed=(mode=="directed"))
    weight <- list(el[,3])
    names(weight) <- weighted
    res <- add_edges(res, edges=t(as.matrix(el[,1:2])), attr=weight)
  } else {
    edges <- unlist(apply(el, 1, function(x) rep(unname(x[1:2]), x[3])))
    res <- graph(n=vc, edges, directed=(mode=="directed"))
  }
  res
}



#' Create graphs from adjacency matrices
#'
#' \code{graph_from_adjacency_matrix} is a flexible function for creating \code{igraph}
#' graphs from adjacency matrices.
#'
#' The order of the vertices are preserved, i.e. the vertex corresponding to
#' the first row will be vertex 0 in the graph, etc.
#'
#' \code{graph_from_adjacency_matrix} operates in two main modes, depending on the
#' \code{weighted} argument.
#'
#' If this argument is \code{NULL} then an unweighted graph is created and an
#' element of the adjacency matrix gives the number of edges to create between
#' the two corresponding vertices.  The details depend on the value of the
#' \code{mode} argument: \describe{ \item{list("directed")}{The graph will be
#' directed and a matrix element gives the number of edges between two
#' vertices.} \item{list("undirected")}{This is exactly the same as \code{max},
#' for convenience. Note that it is \emph{not} checked whether the matrix is
#' symmetric.} \item{list("max")}{An undirected graph will be created and
#' \code{max(A(i,j), A(j,i))} gives the number of edges.}
#' \item{list("upper")}{An undirected graph will be created, only the upper
#' right triangle (including the diagonal) is used for the number of edges.}
#' \item{list("lower")}{An undirected graph will be created, only the lower
#' left triangle (including the diagonal) is used for creating the edges.}
#' \item{list("min")}{undirected graph will be created with \code{min(A(i,j),
#' A(j,i))} edges between vertex \code{i} and \code{j}.} \item{list("plus")}{
#' undirected graph will be created with \code{A(i,j)+A(j,i)} edges between
#' vertex \code{i} and \code{j}.} }
#'
#' If the \code{weighted} argument is not \code{NULL} then the elements of the
#' matrix give the weights of the edges (if they are not zero).  The details
#' depend on the value of the \code{mode} argument: \describe{
#' \item{list("directed")}{The graph will be directed and a matrix element
#' gives the edge weights.} \item{list("undirected")}{First we check that the
#' matrix is symmetric. It is an error if not. Then only the upper triangle is
#' used to create a weighted undirected graph.} \item{list("max")}{An
#' undirected graph will be created and \code{max(A(i,j), A(j,i))} gives the
#' edge weights.} \item{list("upper")}{An undirected graph will be created,
#' only the upper right triangle (including the diagonal) is used (for the edge
#' weights).} \item{list("lower")}{An undirected graph will be created, only
#' the lower left triangle (including the diagonal) is used for creating the
#' edges.} \item{list("min")}{An undirected graph will be created,
#' \code{min(A(i,j), A(j,i))} gives the edge weights.} \item{list("plus")}{An
#' undirected graph will be created, \code{A(i,j)+A(j,i)} gives the edge
#' weights.} }
#'
#' @aliases graph.adjacency
#' @param adjmatrix A square adjacency matrix. From igraph version 0.5.1 this
#' can be a sparse matrix created with the \code{Matrix} package.
#' @param mode Character scalar, specifies how igraph should interpret the
#' supplied matrix. See also the \code{weighted} argument, the interpretation
#' depends on that too. Possible values are: \code{directed},
#' \code{undirected}, \code{upper}, \code{lower}, \code{max}, \code{min},
#' \code{plus}. See details below.
#' @param weighted This argument specifies whether to create a weighted graph
#' from an adjacency matrix. If it is \code{NULL} then an unweighted graph is
#' created and the elements of the adjacency matrix gives the number of edges
#' between the vertices. If it is a character constant then for every non-zero
#' matrix entry an edge is created and the value of the entry is added as an
#' edge attribute named by the \code{weighted} argument. If it is \code{TRUE}
#' then a weighted graph is created and the name of the edge attribute will be
#' \code{weight}. See also details below.
#' @param diag Logical scalar, whether to include the diagonal of the matrix in
#' the calculation. If this is \code{FALSE} then the diagonal is zerod out
#' first.
#' @param add.colnames Character scalar, whether to add the column names as
#' vertex attributes. If it is \sQuote{\code{NULL}} (the default) then, if
#' present, column names are added as vertex attribute \sQuote{name}. If
#' \sQuote{\code{NA}} then they will not be added.  If a character constant,
#' then it gives the name of the vertex attribute to add.
#' @param add.rownames Character scalar, whether to add the row names as vertex
#' attributes. Possible values the same as the previous argument. By default
#' row names are not added. If \sQuote{\code{add.rownames}} and
#' \sQuote{\code{add.colnames}} specify the same vertex attribute, then the
#' former is ignored.
#' @return An igraph graph object.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \link{graph} and \code{\link{graph_from_literal}} for other ways to
#' create graphs.
#' @keywords graphs
#' @examples
#'
#' adjm <- matrix(sample(0:1, 100, replace=TRUE, prob=c(0.9,0.1)), nc=10)
#' g1 <- graph_from_adjacency_matrix( adjm )
#' adjm <- matrix(sample(0:5, 100, replace=TRUE,
#'                       prob=c(0.9,0.02,0.02,0.02,0.02,0.02)), nc=10)
#' g2 <- graph_from_adjacency_matrix(adjm, weighted=TRUE)
#' E(g2)$weight
#'
#' ## various modes for weighted graphs, with some tests
#' nzs <- function(x) sort(x [x!=0])
#' adjm <- matrix(runif(100), 10)
#' adjm[ adjm<0.5 ] <- 0
#' g3 <- graph_from_adjacency_matrix((adjm + t(adjm))/2, weighted=TRUE,
#'                       mode="undirected")
#'
#' g4 <- graph_from_adjacency_matrix(adjm, weighted=TRUE, mode="max")
#' all(nzs(pmax(adjm, t(adjm))[upper.tri(adjm)]) == sort(E(g4)$weight))
#'
#' g5 <- graph_from_adjacency_matrix(adjm, weighted=TRUE, mode="min")
#' all(nzs(pmin(adjm, t(adjm))[upper.tri(adjm)]) == sort(E(g5)$weight))
#'
#' g6 <- graph_from_adjacency_matrix(adjm, weighted=TRUE, mode="upper")
#' all(nzs(adjm[upper.tri(adjm)]) == sort(E(g6)$weight))
#'
#' g7 <- graph_from_adjacency_matrix(adjm, weighted=TRUE, mode="lower")
#' all(nzs(adjm[lower.tri(adjm)]) == sort(E(g7)$weight))
#'
#' g8 <- graph_from_adjacency_matrix(adjm, weighted=TRUE, mode="plus")
#' d2 <- function(x) { diag(x) <- diag(x)/2; x }
#' all(nzs((d2(adjm+t(adjm)))[lower.tri(adjm)]) == sort(E(g8)$weight))
#'
#' g9 <- graph_from_adjacency_matrix(adjm, weighted=TRUE, mode="plus", diag=FALSE)
#' d0 <- function(x) { diag(x) <- 0 }
#' all(nzs((d0(adjm+t(adjm)))[lower.tri(adjm)]) == sort(E(g9)$weight))
#'
#' ## row/column names
#' rownames(adjm) <- sample(letters, nrow(adjm))
#' colnames(adjm) <- seq(ncol(adjm))
#' g10 <- graph_from_adjacency_matrix(adjm, weighted=TRUE, add.rownames="code")
#' summary(g10)
#'
graph_from_adjacency_matrix <- function(adjmatrix, mode=c("directed", "undirected", "max",
                                         "min", "upper", "lower", "plus"),
                            weighted=NULL, diag=TRUE,
                            add.colnames=NULL, add.rownames=NA) {

  if (inherits(adjmatrix, "Matrix")) {
    res <- graph.adjacency.sparse(adjmatrix, mode=mode, weighted=weighted, diag=diag)
  } else {
    res <- graph.adjacency.dense(adjmatrix, mode=mode, weighted=weighted, diag=diag)
  }

  ## Add columns and row names as attributes
  if (is.null(add.colnames)) {
    if (!is.null(colnames(adjmatrix))) {
      add.colnames <- "name"
    } else {
      add.colnames <- NA
    }
  } else if (!is.na(add.colnames)) {
    if (is.null(colnames(adjmatrix))) {
      warning("No column names to add")
      add.colnames <- NA
    }
  }

  if (is.null(add.rownames)) {
    if (!is.null(rownames(adjmatrix))) {
      add.rownames <- "name"
    } else {
      add.colnames <- NA
    }
  } else if (!is.na(add.rownames)) {
    if (is.null(rownames(adjmatrix))) {
      warning("No row names to add")
      add.rownames <- NA
    }
  }

  if (!is.na(add.rownames) && !is.na(add.colnames) &&
      add.rownames == add.colnames ) {
    warning("Same attribute for columns and rows, row names are ignored")
    add.rownames <- NA
  }

  if (!is.na(add.colnames)) {
    res <- set_vertex_attr(res, add.colnames, value=colnames(adjmatrix))
  }
  if (!is.na(add.rownames)) {
    res <- set_vertex_attr(res, add.rownames, value=rownames(adjmatrix))
  }

  res
}

#' @rdname graph_from_adjacency_matrix
#' @param ... Passed to \code{graph_from_adjacency_matrix}.
#' @export

from_adjacency <- function(...) constructor_spec(graph_from_adjacency_matrix, ...)
