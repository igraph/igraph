#   IGraph R package
#   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#' Create an igraph graph from a list of edges
#'
#' @param edges Numeric vector defining the edges, the first edge points
#'   from the first element to the second, the second edge from the third
#'   to the fourth, etc.
#' @param n The number of vertices in the graph. This parameter is ignored
#'   if there is a bigger vertex id in \code{edges}. This means that for
#'   this function it is safe to supply zero here if the vertex with the
#'   largest id is not an isolate.
#' @param directed Whether to create a directed graph.
#' @return An igraph graph.
#'
#' @family determimistic constructors
#' @export
#' @examples
#' graph(c(1, 2, 2, 3, 3, 4, 5, 6), directed = FALSE)

graph <- function(edges, n=max(edges), directed=TRUE ) {
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_create", as.numeric(edges)-1, as.numeric(n),
        as.logical(directed),
        PACKAGE="igraph")
}

#' A graph with no edges
#'
#' @aliases graph.empty
#' @concept Empty graph.
#' @param n Number of vertices.
#' @param directed Whether to create a directed graph.
#' @return An igraph graph.
#'
#' @family determimistic constructors
#' @export
#' @examples
#' empty_graph(n = 10)
#' empty_graph(n = 5, directed = FALSE)

empty_graph <- function(n=0, directed=TRUE) {
  # Argument checks
  n <- as.integer(n)
  directed <- as.logical(directed)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_empty", n, directed,
        PACKAGE="igraph")

  res
}

#' Creating (small) graphs via a simple interface
#'
#' This function is useful if you want to create a small (named) graph
#' quickly, it works for both directed and undirected graphs.
#'
#' @details
#' \code{graph_from_formula} is very handy for creating small graphs quickly.
#' You need to supply one or more R expressions giving the structure of
#' the graph. The expressions consist of vertex names and edge
#' operators. An edge operator is a sequence of \sQuote{\code{-}} and
#' \sQuote{\code{+}} characters, the former is for the edges and the
#' latter is used for arrow heads. The edges can be arbitrarily long,
#' ie. you may use as many \sQuote{\code{-}} characters to \dQuote{draw}
#' them as you like.
#'
#' If all edge operators consist of only \sQuote{\code{-}} characters
#' then the graph will be undirected, whereas a single \sQuote{\code{+}}
#' character implies a directed graph.
#'
#' Let us see some simple examples. Without arguments the function
#' creates an empty graph:
#' \preformatted{  graph_from_formula()
#' }
#'
#' A simple undirected graph with two vertices called \sQuote{A} and
#' \sQuote{B} and one edge only:
#' \preformatted{  graph_from_formula(A-B)
#' }
#'
#' Remember that the length of the edges does not matter, so we could
#' have written the following, this creates the same graph:
#' \preformatted{  graph_from_formula( A-----B )
#' }
#'
#' If you have many disconnected components in the graph, separate them
#' with commas. You can also give isolate vertices.
#' \preformatted{  graph_from_formula( A--B, C--D, E--F, G--H, I, J, K )
#' }
#'
#' The \sQuote{\code{:}} operator can be used to define vertex sets. If
#' an edge operator connects two vertex sets then every vertex from the
#' first set will be connected to every vertex in the second set. The
#' following form creates a full graph, including loop edges:
#' \preformatted{  graph_from_formula( A:B:C:D -- A:B:C:D )
#' }
#'
#' In directed graphs, edges will be created only if the edge operator
#' includes a arrow head (\sQuote{+}) \emph{at the end} of the edge:
#' \preformatted{  graph_from_formula( A -+ B -+ C )
#'   graph_from_formula( A +- B -+ C )
#'   graph_from_formula( A +- B -- C )
#' }
#' Thus in the third example no edge is created between vertices \code{B}
#' and \code{C}.
#'
#' Mutual edges can be also created with a simple edge operator:
#' \preformatted{  graph_from_formula( A +-+ B +---+ C ++ D + E)
#' }
#' Note again that the length of the edge operators is arbitrary,
#' \sQuote{\code{+}}, \sQuote{\code{++}} and \sQuote{\code{+-----+}} have
#' exactly the same meaning.
#'
#' If the vertex names include spaces or other special characters then
#' you need to quote them:
#' \preformatted{  graph_from_formula( "this is" +- "a silly" -+ "graph here" )
#' }
#' You can include any character in the vertex names this way, even
#' \sQuote{+} and \sQuote{-} characters.
#'
#' See more examples below.
#'
#' @aliases graph.formula
#' @param ... The formulae giving the structure of the graph, see
#'   details below.
#' @param simplify Logical scalar, whether to call \code{\link{simplify}}
#'   on the created graph. By default the graph is simplified, loop and
#'   multiple edges are removed.
#' @return An igraph graph
#'
#' @family determimistic constructors
#' @export
#' @examples
#' # A simple undirected graph
#' g <- graph_from_formula( Alice-Bob-Cecil-Alice, Daniel-Cecil-Eugene,
#'                      Cecil-Gordon )
#' g
#'
#' # Another undirected graph, ":" notation
#' g2 <- graph_from_formula( Alice-Bob:Cecil:Daniel, Cecil:Daniel-Eugene:Gordon )
#' g2
#'
#' # A directed graph
#' g3 <- graph_from_formula( Alice +-+ Bob --+ Cecil +-- Daniel,
#'                      Eugene --+ Gordon:Helen )
#' g3
#'
#' # A graph with isolate vertices
#' g4 <- graph_from_formula( Alice -- Bob -- Daniel, Cecil:Gordon, Helen )
#' g4
#' V(g4)$name
#'
#' # "Arrows" can be arbitrarily long
#' g5 <- graph_from_formula( Alice +---------+ Bob )
#' g5
#'
#' # Special vertex names
#' g6 <- graph_from_formula( "+" -- "-", "*" -- "/", "%%" -- "%/%" )
#' g6
#'

graph_from_formula <- function(..., simplify=TRUE) {
  mf <- as.list(match.call())[-1]

  ## In case 'simplify' is given
  if ('simplify' %in% names(mf)) {
    w <- which(names(mf)=='simplify')
    if (length(w) > 1) { stop("'simplify' specified multiple times") }
    mf <- mf[-w]
  }

  ## Operators first
  f <- function(x) {
    if (is.call(x)) {
      return (list(as.character(x[[1]]), lapply(x[-1], f)))
    } else {
      return (NULL)
    }
  }
  ops <- unlist(lapply(mf, f))
  if (all(ops %in% c("-", ":"))) {
    directed <- FALSE
  } else if (all(ops %in% c("-", "+", ":"))) {
    directed <- TRUE
  } else {
    stop("Invalid operator in formula")
  }

  f <- function(x) {
    if (is.call(x)) {
      if (length(x)==3) {
        return( list(f(x[[2]]), op=as.character(x[[1]]), f(x[[3]])) )
      } else {
        return( list(op=as.character(x[[1]]), f(x[[2]])) )
      }
    } else {
      return( c(sym=as.character(x)) )
    }
  }

  ret <- lapply(mf, function(x) unlist(f(x)))

  v <- unique(unlist(lapply(ret, function(x) { x[ names(x)=="sym" ] })))

  ## Merge symbols for ":"
  ret <- lapply(ret, function(x) {
    res <- list()
    for (i in seq(along=x)) {
      if (x[i]==":" && names(x)[i]=="op") {
        ## SKIP
      } else if (i>1 && x[i-1]==":" && names(x)[i-1]=="op") {
        res[[length(res)]] <- c(res[[length(res)]], unname(x[i]))
      } else {
        res <- c(res, x[i])
      }
    }
    res
  })

  ## Ok, create the edges
  edges <- numeric()
  for (i in seq(along=ret)) {
    prev.sym <- character()
    lhead <- rhead <- character()
    for (j in seq(along=ret[[i]])) {
      act <- ret[[i]][[j]]
      if (names(ret[[i]])[j]=="op") {
        if (length(lhead)==0) {
          lhead <- rhead <- act
        } else {
          rhead <- act
        }
      } else if (names(ret[[i]])[j]=="sym") {
        for (ps in prev.sym) {
          for (ps2 in act) {
            if (lhead=="+") {
              edges <- c(edges, unname(c(ps2, ps)))
            }
            if (!directed || rhead=="+") {
              edges <- c(edges, unname(c(ps, ps2)))
            }
          }
        }
        lhead <- rhead <- character()
        prev.sym <- act
      }
    }
  }

  ids <- seq(along=v)
  names(ids) <- v
  res <- graph( unname(ids[edges]), n=length(v), directed=directed)
  if (simplify) res <- simplify(res)
  res <- set_vertex_attr(res, "name", value=v)
  res
}

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
    res <- empty_graph(n=vc, directed=(mode=="directed"))
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
#' @seealso \link{graph} and \code{\link{graph_from_formula}} for other ways to
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


#' Create a star graph, a tree with n vertices and n - 1 leaves
#'
#' \code{star} creates a star graph, in this every single vertex is
#' connected to the center vertex and nobody else.
#'
#' @aliases graph.star
#' @concept Star graph
#' @param n Number of vertices.
#' @param mode It defines the direction of the
#'   edges, \code{in}: the edges point \emph{to} the center, \code{out}:
#'   the edges point \emph{from} the center, \code{mutual}: a directed
#'   star is created with mutual edges, \code{undirected}: the edges
#'   are undirected.
#' @param center ID of the center vertex.
#' @return An igraph graph.
#'
#' @family determimistic constructors
#' @export
#' @examples
#' star(10, mode = "out")
#' star(5, mode = "undirected")

star <- function(n, mode=c("in", "out", "mutual", "undirected"),
                   center=1 ) {

  mode <- igraph.match.arg(mode)
  mode1 <- switch(mode, "out"=0, "in"=1, "undirected"=2, "mutual"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_star", as.numeric(n), as.numeric(mode1),
               as.numeric(center)-1,
               PACKAGE="igraph")
  if (getIgraphOpt("add.params")) {
    res$name <- switch(mode, "in"="In-star", "out"="Out-star", "Star")
    res$mode <- mode
    res$center <- center
  }
  res
}

#' Create a full graph
#'
#' @aliases graph.full
#' @concept Full graph
#' @param n Number of vertices.
#' @param directed Whether to create a directed graph.
#' @param loops Whether to add self-loops to the graph.
#' @return An igraph graph
#'
#' @family determimistic constructors
#' @export
#' @examples
#' full_graph(5)
#' str(full_graph(4, directed = TRUE))

full_graph <- function(n, directed=FALSE, loops=FALSE) {
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_full", as.numeric(n), as.logical(directed),
               as.logical(loops),
               PACKAGE="igraph")
  if (getIgraphOpt("add.params")) {
    res$name <- "Full graph"
    res$loops <- loops
  }
  res
}

###################################################################
# Lattices, every kind
###################################################################

#' Create a lattice graph
#'
#' \code{lattice} is a flexible function, it can create lattices of
#' arbitrary dimensions, periodic or unperiodic ones. It has two
#' forms. In the first form you only supply \code{dimvector}, but not
#' \code{length} and \code{dim}. In the second form you omit
#' \code{dimvector} and supply \code{length} and \code{dim}.
#'
#' @aliases graph.lattice
#' @concept Lattice
#' @param dimvector A vector giving the size of the lattice in each
#'   dimension.
#' @param length Integer constant, for regular lattices, the size of the
#'   lattice in each dimension.
#' @param dim Integer constant, the dimension of the lattice.
#' @param nei The distance within which (inclusive) the neighbors on the
#'   lattice will be connected. This parameter is not used right now.
#' @param directed Whether to create a directed lattice.
#' @param mutual Logical, if \code{TRUE} directed lattices will be
#'   mutually connected.
#' @param circular Logical, if \code{TRUE} the lattice or ring will be
#'   circular.
#' @param ... Currently ignored.
#' @return An igraph graph.
#'
#' @family determimistic constructors
#' @export
#' @examples
#' lattice(c(5, 5, 5))
#' lattice(length = 5, dim = 3)

lattice <- function(dimvector = NULL, length = NULL, dim = NULL,
                          nei = 1, directed = FALSE, mutual = FALSE,
                          circular=FALSE, ...) {

##   # Check
##   if (is.null(dimvector) && (is.null(length) || is.null(dim))) {
##     stop("Either `length' and `dim' or 'dimvector' must be set. See docs.")
##   }
##   if (!is.null(length) && length < 1) {
##     stop("Invalid `length' argument, should be at least one")
##   }
##   if (!is.null(length) && dim < 1) {
##     stop("Invalid `dim' argument, should be at least one")
##   }
##   if (!is.null(length) && any(dimvector < 1)) {
##     stop("Invalid `dimvector', has negative or smaller than one elements")
##   }
##   if (mutual && !directed) {
##     warning("`mutual' specified for undirected graph, proceeding with multiplex edges...")
##   }
##   if (nei < 1) {
##     stop("`nei' should be at least one")
##   }

##   if (!is.null(length)) {
##     length <- as.numeric(length)
##     dim <- as.numeric(dim)
##     dimvector <- rep(length, times=dim)
##   } else {
##     dimvector <- as.numeric(dimvector)
##   }
##   nei <- as.numeric(nei)

##   n <- prod(dimvector)
##   res <- empty_graph(n=n, directed=directed, ...)
##   res <- add_edges(res, .Call("REST_create_lattice", dimvector, n,
##                               circular, mutual, PACKAGE="igraph"))

##   # Connect also to local neighborhood
##   if (nei >= 2) {
##     neighbors <- lapply(1:length(res), function(a) get.neighborhood(res, a))
##     res <- add_edges(res, .Call("REST_connect_neighborhood", neighbors, nei,
##                                 mutual, PACKAGE="igraph"))
##   }

##   res

  if (is.null(dimvector)) {
    dimvector <- rep(length, dim)
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_lattice", as.numeric(dimvector), as.numeric(nei),
               as.logical(directed), as.logical(mutual),
               as.logical(circular),
               PACKAGE="igraph")
  if (getIgraphOpt("add.params")) {
    res$name <- "Lattice graph"
    res$dimvector <- dimvector
    res$nei <- nei
    res$mutual <- mutual
    res$circular <- circular
  }
  res
}

#' Create a ring graph
#'
#' A ring is a one-dimensional lattice and this function is a special case
#' of \code{\link{lattice}}.
#'
#' @aliases ring graph.ring
#' @param n Number of vertices.
#' @param directed Whether the graph is directed.
#' @param mutual Whether directed edges are mutual. It is ignored in
#'   undirected graphs.
#' @param circular Whether to create a circular ring. A non-circular
#'   ring is essentially a \dQuote{line}: a tree where every non-leaf
#'   vertex has one child.
#' @return An igraph graph.
#'
#' @family determimistic constructors
#' @export
#' @examples
#' str(ring(10))
#' str(ring(10, directed = TRUE, mutual = TRUE))

ring <- function(n, directed=FALSE, mutual=FALSE, circular=TRUE) {
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_ring", as.numeric(n), as.logical(directed),
               as.logical(mutual), as.logical(circular),
               PACKAGE="igraph")
  if (getIgraphOpt("add.params")) {
    res$name <- "Ring graph"
    res$mutual <- mutual
    res$circular <- circular
  }
  res
}


###################################################################
# Trees, regular
###################################################################

#' Create tree graphs
#'
#' Create a regular tree graph.
#'
#' @aliases graph.tree
#' @concept Trees.
#' @param n Number of vertices.
#' @param children Integer scalar, the number of children of a vertex
#'   (except for leafs)
#' @param mode Defines the direction of the
#'   edges. \code{out} indicates that the edges point from the parent to
#'   the children, \code{in} indicates that they point from the children
#'   to their parents, while \code{undirected} creates an undirected
#'   graph.
#' @return An igraph graph
#'
#' @family determimistic constructors
#' @export
#' @examples
#' tree(10, 2)
#' tree(10, 3, mode = "undirected")

tree <- function(n, children=2, mode=c("out", "in", "undirected")) {

  mode <- igraph.match.arg(mode)
  mode1 <- switch(mode, "out"=0, "in"=1, "undirected"=2);

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_tree", as.numeric(n), as.numeric(children),
               as.numeric(mode1),
               PACKAGE="igraph")
  if (getIgraphOpt("add.params")) {
    res$name <- "Tree"
    res$children <- children
    res$mode <- mode
  }
  res
}

###################################################################
# The graph atlas
###################################################################

#' Create a graph from the Graph Atlas
#'
#' \code{graph_atlas} creates graphs from the book
#' \sQuote{An Atlas of Graphs} by
#' Roland C. Read and Robin J. Wilson. The atlas contains all undirected
#' graphs with up to seven vertices, numbered from 0 up to 1252. The
#' graphs are listed:
#' \enumerate{
#'    \item in increasing order of number of nodes;
#'    \item for a fixed number of nodes, in increasing order of the number
#'      of edges;
#'    \item for fixed numbers of nodes and edges, in increasing order of
#'      the degree sequence, for example 111223 < 112222;
#'    \item for fixed degree sequence, in increasing number of
#'      automorphisms.
#' }
#'
#' @aliases graph.atlas
#' @concept Graph Atlas.
#' @param n The id of the graph to create.
#' @return An igraph graph.
#'
#' @family determimistic constructors
#' @export
#' @examples
#' ## Some randomly picked graphs from the atlas
#' graph_atlas(sample(0:1252, 1))
#' graph_atlas(sample(0:1252, 1))


graph_atlas <- function(n) {

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_atlas", as.numeric(n),
               PACKAGE="igraph")
  if (getIgraphOpt("add.params")) {
    res$name <- sprintf("Graph from the Atlas #%i", n)
    res$n <- n
  }
  res
}

###################################################################
# Create a graph from a data frame
###################################################################



#' Creating igraph graphs from data frames or vice-versa
#' 
#' This function creates an igraph graph from one or two data frames containing
#' the (symbolic) edge list and edge/vertex attributes.
#' 
#' \code{graph_from_data_frame} creates igraph graphs from one or two data frames.
#' It has two modes of operatation, depending whether the \code{vertices}
#' argument is \code{NULL} or not.
#' 
#' If \code{vertices} is \code{NULL}, then the first two columns of \code{d}
#' are used as a symbolic edge list and additional columns as edge attributes.
#' The names of the attributes are taken from the names of the columns.
#' 
#' If \code{vertices} is not \code{NULL}, then it must be a data frame giving
#' vertex metadata. The first column of \code{vertices} is assumed to contain
#' symbolic vertex names, this will be added to the graphs as the
#' \sQuote{\code{name}} vertex attribute. Other columns will be added as
#' additional vertex attributes. If \code{vertices} is not \code{NULL} then the
#' symbolic edge list given in \code{d} is checked to contain only vertex names
#' listed in \code{vertices}.
#' 
#' Typically, the data frames are exported from some speadsheat software like
#' Excel and are imported into R via \code{\link{read.table}},
#' \code{\link{read.delim}} or \code{\link{read.csv}}.
#' 
#' \code{as_data_frame} converts the igraph graph into one or more data
#' frames, depending on the \code{what} argument.
#' 
#' If the \code{what} argument is \code{edges} (the default), then the edges of
#' the graph and also the edge attributes are returned. The edges will be in
#' the first two columns, named \code{from} and \code{to}. (This also denotes
#' edge direction for directed graphs.)  For named graphs, the vertex names
#' will be included in these columns, for other graphs, the numeric vertex ids.
#' The edge attributes will be in the other columns. It is not a good idea to
#' have an edge attribute named \code{from} or \code{to}, because then the
#' column named in the data frame will not be unique. The edges are listed in
#' the order of their numeric ids.
#' 
#' If the \code{what} argument is \code{vertices}, then vertex attributes are
#' returned. Vertices are listed in the order of their numeric vertex ids.
#' 
#' If the \code{what} argument is \code{both}, then both vertex and edge data
#' is returned, in a list with named entries \code{vertices} and \code{edges}.
#' 
#' @aliases graph_from_data_frame graph.data.frame as_data_frame get.data.frame
#' @param d A data frame containing a symbolic edge list in the first two
#' columns. Additional columns are considered as edge attributes.  Since
#' version 0.7 this argument is coerced to a data frame with
#' \code{as.data.frame}.
#' @param directed Logical scalar, whether or not to create a directed graph.
#' @param vertices A data frame with vertex metadata, or \code{NULL}. See
#' details below. Since version 0.7 this argument is coerced to a data frame
#' with \code{as.data.frame}, if not \code{NULL}.
#' @return An igraph graph object for \code{graph_from_data_frame}, and either a
#' data frame or a list of two data frames named \code{edges} and
#' \code{vertices} for \code{as.data.frame}.
#' @note For \code{graph_from_data_frame} \code{NA} elements in the first two
#' columns \sQuote{d} are replaced by the string \dQuote{NA} before creating
#' the graph. This means that all \code{NA}s will correspond to a single
#' vertex.
#' 
#' \code{NA} elements in the first column of \sQuote{vertices} are also
#' replaced by the string \dQuote{NA}, but the rest of \sQuote{vertices} is not
#' touched. In other words, vertex names (=the first column) cannot be
#' \code{NA}, but other vertex attributes can.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{graph_from_formula}}
#' for another way to create graphs, \code{\link{read.table}} to read in tables
#' from files.
#' @keywords graphs
#' @examples
#' 
#' ## A simple example with a couple of actors
#' ## The typical case is that these tables are read in from files....
#' actors <- data.frame(name=c("Alice", "Bob", "Cecil", "David",
#'                             "Esmeralda"),
#'                      age=c(48,33,45,34,21),
#'                      gender=c("F","M","F","M","F"))
#' relations <- data.frame(from=c("Bob", "Cecil", "Cecil", "David",
#'                                "David", "Esmeralda"),
#'                         to=c("Alice", "Bob", "Alice", "Alice", "Bob", "Alice"),
#'                         same.dept=c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),
#'                         friendship=c(4,5,5,2,1,1), advice=c(4,5,5,4,2,3))
#' g <- graph_from_data_frame(relations, directed=TRUE, vertices=actors)
#' print(g, e=TRUE, v=TRUE)
#' 
#' ## The opposite operation
#' as_data_frame(g, what="vertices")
#' as_data_frame(g, what="edges")
#' 
graph_from_data_frame <- function(d, directed=TRUE, vertices=NULL) {

  d <- as.data.frame(d)
  if (!is.null(vertices)) { vertices <- as.data.frame(vertices) }

  if (ncol(d) < 2) {
    stop("the data frame should contain at least two columns")
  }

  ## Handle if some elements are 'NA'
  if (any(is.na(d[,1:2]))) {
    warning("In `d' `NA' elements were replaced with string \"NA\"")
    d[,1:2][ is.na(d[,1:2]) ] <- 'NA'
  }
  if (!is.null(vertices) && any(is.na(vertices[,1]))) {
    warning("In `vertices[,1]' `NA' elements were replaced with string \"NA\"")
    vertices[,1][is.na(vertices[,1])] <- 'NA'
  }

  names <- unique( c(as.character(d[,1]), as.character(d[,2])) )
  if (!is.null(vertices)) {
    names2 <- names
    vertices <- as.data.frame(vertices)
    if (ncol(vertices) < 1) {
      stop("Vertex data frame contains no rows")
    }
    names <- as.character(vertices[,1])
    if (any(duplicated(names))) {
      stop("Duplicate vertex names")
    }
    if (any(! names2 %in% names)) {
      stop("Some vertex names in edge list are not listed in vertex data frame")
    }
  }

  # create graph
  g <- empty_graph(n=0, directed=directed)

  # vertex attributes
  attrs <- list(name=names)
  if (!is.null(vertices)) {
    if (ncol(vertices) > 1) {
      for (i in 2:ncol(vertices)) {
        newval <- vertices[,i]
        if (class(newval) == "factor") {
          newval <- as.character(newval)
        }
        attrs[[ names(vertices)[i] ]] <- newval
      }
    }
  }

  # add vertices
  g <- add_vertices(g, length(names), attr=attrs)

  # create edge list
  from <- as.character(d[,1])
  to <- as.character(d[,2])
  edges <- rbind(match(from, names), match(to,names))

  # edge attributes
  attrs <- list()
  if (ncol(d) > 2) {
    for (i in 3:ncol(d)) {
      newval <- d[,i]
      if (class(newval) == "factor") {
        newval <- as.character(newval)
      }
      attrs[[ names(d)[i] ]] <- newval
    }
  }

  # add the edges
  g <- add_edges(g, edges, attr=attrs)
  g
}

#' Create a graph from an edge list matrix
#'
#' \code{graph_from_edgelist} creates a graph from an edge list. Its argument
#' is a two-column matrix, each row defines one edge. If it is
#' a numeric matrix then its elements are interpreted as vertex ids. If
#' it is a character matrix then it is interpreted as symbolic vertex
#' names and a vertex id will be assigned to each name, and also a
#' \code{name} vertex attribute will be added.
#'
#' @aliases graph.edgelist
#' @concept Edge list
#' @param el The edge list, a two column matrix, character or numeric.
#' @param directed Whether to create a directed graph.
#' @return An igraph graph.
#'
#' @family determimistic constructors
#' @export
#' @examples
#' el <- matrix( c("foo", "bar", "bar", "foobar"), nc = 2, byrow = TRUE)
#' graph_from_edgelist(el)
#'
#' # Create a ring by hand
#' graph_from_edgelist(cbind(1:10, c(2:10, 1)))

graph_from_edgelist <- function(el, directed=TRUE) {

  if (!is.matrix(el) || ncol(el) != 2) {
    stop("graph_from_edgelist expects a matrix with two columns")
  }

  if (nrow(el) == 0) {
    res <- empty_graph(directed=directed)
  } else {
    if (is.character(el)) {
      ## symbolic edge list
      names <- unique(as.character(t(el)))
      ids <- seq(names)
      names(ids) <- names
      res <- graph( unname(ids[t(el)]), directed=directed)
      rm(ids)
      V(res)$name <- names
    } else {
      ## normal edge list
      res <- graph( t(el), directed=directed )
    }
  }

  res
}

#' Create an extended chordal ring graph
#'
#' \code{chordal_ring} creates an extended chordal ring.
#' An extended chordal ring is regular graph, each node has the same
#' degree. It can be obtained from a simple ring by adding some extra
#' edges specified by a matrix. Let p denote the number of columns in
#' the \sQuote{\code{W}} matrix. The extra edges of vertex \code{i}
#' are added according to column \code{i mod p} in
#' \sQuote{\code{W}}. The number of extra edges is the number
#' of rows in \sQuote{\code{W}}: for each row \code{j} an edge
#' \code{i->i+w[ij]} is added if \code{i+w[ij]} is less than the number
#' of total nodes. See also Kotsis, G: Interconnection Topologies for
#' Parallel Processing Systems, PARS Mitteilungen 11, 1-6, 1993.
#'
#' @aliases graph.extended.chordal.ring
#' @param n The number of vertices.
#' @param w A matrix which specifies the extended chordal ring. See
#'   details below.
#' @return An igraph graph.
#'
#' @family determimistic constructors
#' @export
#' @examples
#' chord <- chordal_ring(15,
#'     matrix(c(3, 12, 4, 7, 8, 11), nr = 2))

chordal_ring <- function(n, w) {

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_extended_chordal_ring", as.numeric(n),
               as.matrix(w),
               PACKAGE="igraph")
  if (getIgraphOpt("add.params")) {
    res$name <- "Extended chordal ring"
    res$w <- w
  }
  res
}



#' Line graph of a graph
#' 
#' This function calculates the line graph of another graph.
#' 
#' The line graph \code{L(G)} of a \code{G} undirected graph is defined as
#' follows. \code{L(G)} has one vertex for each edge in \code{G} and two
#' vertices in \code{L(G)} are connected by an edge if their corresponding
#' edges share an end point.
#' 
#' The line graph \code{L(G)} of a \code{G} directed graph is slightly
#' different, \code{L(G)} has one vertex for each edge in \code{G} and two
#' vertices in \code{L(G)} are connected by a directed edge if the target of
#' the first vertex's corresponding edge is the same as the source of the
#' second vertex's corresponding edge.
#'
#' @aliases line.graph
#' @param graph The input graph, it can be directed or undirected.
#' @return A new graph object.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}, the first version of
#' the C code was written by Vincent Matossian.
#' @keywords graphs
#' @examples
#' 
#' # generate the first De-Bruijn graphs
#' g <- full_graph(2, directed=TRUE, loops=TRUE)
#' line_graph(g)
#' line_graph(line_graph(g))
#' line_graph(line_graph(line_graph(g)))
#' 
line_graph <- function(graph) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_line_graph", graph,
               PACKAGE="igraph")
  if (getIgraphOpt("add.params")) {
    res$name <- "Line graph"
  }
  res
}




#' De Bruijn graphs.
#' 
#' De Bruijn graphs are labeled graphs representing the overlap of strings.
#' 
#' A de Bruijn graph represents relationships between strings. An alphabet of
#' \code{m} letters are used and strings of length \code{n} are considered.  A
#' vertex corresponds to every possible string and there is a directed edge
#' from vertex \code{v} to vertex \code{w} if the string of \code{v} can be
#' transformed into the string of \code{w} by removing its first letter and
#' appending a letter to it.
#' 
#' Please note that the graph will have \code{m} to the power \code{n} vertices
#' and even more edges, so probably you don't want to supply too big numbers
#' for \code{m} and \code{n}.
#' 
#' De Bruijn graphs have some interesting properties, please see another
#' source, eg. Wikipedia for details.
#'
#' @aliases graph.de.bruijn
#' @param m Integer scalar, the size of the alphabet. See details below.
#' @param n Integer scalar, the length of the labels. See details below.
#' @return A graph object.
#' @author Gabor Csardi <csardi.gabor@@gmail.com>
#' @seealso \code{\link{kautz_graph}}, \code{\link{line_graph}}
#' @keywords graphs
#' @examples
#' 
#' # de Bruijn graphs can be created recursively by line graphs as well 
#' g <- de_bruijn_graph(2,1)
#' de_bruijn_graph(2,2)
#' line_graph(g)
#' 
de_bruijn_graph <- function(m, n) {

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_de_bruijn", as.numeric(m), as.numeric(n),
               PACKAGE="igraph")
  if (getIgraphOpt("add.params")) {
    res$name <- sprintf("De-Bruijn graph %i-%i", m, n)
    res$m <- m
    res$n <- n
  }
  res
}



#' Kautz graphs
#' 
#' Kautz graphs are labeled graphs representing the overlap of strings.
#' 
#' A Kautz graph is a labeled graph, vertices are labeled by strings of length
#' \code{n+1} above an alphabet with \code{m+1} letters, with the restriction
#' that every two consecutive letters in the string must be different. There is
#' a directed edge from a vertex \code{v} to another vertex \code{w} if it is
#' possible to transform the string of \code{v} into the string of \code{w} by
#' removing the first letter and appending a letter to it.
#' 
#' Kautz graphs have some interesting properties, see eg. Wikipedia for
#' details.
#'
#' @aliases graph.kautz
#' @param m Integer scalar, the size of the alphabet. See details below.
#' @param n Integer scalar, the length of the labels. See details below.
#' @return A graph object.
#' @author Gabor Csardi <csardi.gabor@@gmail.com>, the first version in R was
#' written by Vincent Matossian.
#' @seealso \code{\link{de_bruijn_graph}}, \code{\link{line_graph}}
#' @keywords graphs
#' @examples
#' 
#' line_graph(kautz_graph(2,1))
#' kautz_graph(2,2)
#' 
kautz_graph <- function(m, n) {

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_kautz", as.numeric(m), as.numeric(n),
               PACKAGE="igraph")
  if (getIgraphOpt("add.params")) {
    res$name <- sprintf("Kautz graph %i-%i", m, n)
    res$m <- m
    res$n <- n
  }
  res
}

#' Create a full bipartite graph
#' 
#' Bipartite graphs are also called two-mode by some. This function creates a
#' bipartite graph in which every possible edge is present.
#' 
#' Bipartite graphs have a \sQuote{\code{type}} vertex attribute in igraph,
#' this is boolean and \code{FALSE} for the vertices of the first kind and
#' \code{TRUE} for vertices of the second kind.
#'
#' @aliases graph.full.bipartite
#' @param n1 The number of vertices of the first kind.
#' @param n2 The number of vertices of the second kind.
#' @param directed Logical scalar, whether the graphs is directed.
#' @param mode Scalar giving the kind of edges to create for directed graphs.
#' If this is \sQuote{\code{out}} then all vertices of the first kind are
#' connected to the others; \sQuote{\code{in}} specifies the opposite
#' direction; \sQuote{\code{all}} creates mutual edges. This argument is
#' ignored for undirected graphs.x
#' @return An igraph graph, with the \sQuote{\code{type}} vertex attribute set.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{full_graph}} for creating one-mode full graphs
#' @keywords graphs
#' @examples
#' 
#' g <- full_bipartite_graph(2, 3)
#' g2 <- full_bipartite_graph(2, 3, dir=TRUE)
#' g3 <- full_bipartite_graph(2, 3, dir=TRUE, mode="in")
#' g4 <- full_bipartite_graph(2, 3, dir=TRUE, mode="all")
#' 
full_bipartite_graph <- function(n1, n2, directed=FALSE,
                       mode=c("all", "out", "in")) {

  n1 <- as.integer(n1)
  n2 <- as.integer(n2)
  directed <- as.logical(directed)
  mode1 <- switch(igraph.match.arg(mode), "out"=1, "in"=2, "all"=3, "total"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_full_bipartite", n1, n2, as.logical(directed), mode1,
               PACKAGE="igraph")
  if (getIgraphOpt("add.params")) {
    res$graph$name <- "Full bipartite graph"
    res$n1 <- n1
    res$n2 <- n2
    res$mode <- mode
  }
  set_vertex_attr(res$graph, "type", value=res$types)
}



#' Create a bipartite graph
#' 
#' A bipartite graph has two kinds of vertices and connections are only allowed
#' between different kinds.
#' 
#' Bipartite graphs have a \code{type} vertex attribute in igraph, this is
#' boolean and \code{FALSE} for the vertices of the first kind and \code{TRUE}
#' for vertices of the second kind.
#' 
#' \code{bipartite_graph} basically does three things. First it checks tha
#' \code{edges} vector against the vertex \code{types}. Then it creates a graph
#' using the \code{edges} vector and finally it adds the \code{types} vector as
#' a vertex attribute called \code{type}.
#' 
#' \code{is_bipartite} checks whether the graph is bipartite or not. It just
#' checks whether the graph has a vertex attribute called \code{type}.
#' 
#' @aliases bipartite_graph graph.bipartite is.bipartite is_bipartite
#' @param types A vector giving the vertex types. It will be coerced into
#' boolean. The length of the vector gives the number of vertices in the graph.
#' @param edges A vector giving the edges of the graph, the same way as for the
#' regular \code{\link{graph}} function. It is checked that the edges indeed
#' connect vertices of different kind, accoding to the supplied \code{types}
#' vector.
#' @param directed Whether to create a directed graph, boolean constant. Note
#' that by default undirected graphs are created, as this is more common for
#' bipartite graphs.
#' @param graph The input graph.
#' @return \code{bipartite_graph} returns a bipartite igraph graph. In other
#' words, an igraph graph that has a vertex attribute named \code{type}.
#' 
#' \code{is_bipartite} returns a logical scalar.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{graph}} to create one-mode networks
#' @keywords graphs
#' @examples
#' 
#' g <- bipartite_graph( rep(0:1,length=10), c(1:10))
#' print(g, v=TRUE)
#' 
bipartite_graph <- function(types, edges, directed=FALSE) {

  types <- as.logical(types)
  edges <- as.numeric(edges)-1
  directed <- as.logical(directed)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_create_bipartite", types, edges, directed,
               PACKAGE="igraph")
  set_vertex_attr(res, "type", value=types)
}

graph.incidence.sparse <- function(incidence, directed, mode, multiple,
                                   weighted) {
  n1 <- nrow(incidence)
  n2 <- ncol(incidence)
  el <- Matrix::summary(incidence)
  ## el <- summary(incidence)
  el[,2] <- el[,2] + n1

  if (!is.null(weighted)) {

    if (is.logical(weighted) && weighted) {
      weighted <- "weight"
    }
    if (!is.character(weighted)) {
      stop("invalid value supplied for `weighted' argument, please see docs.")
    }

    if (!directed || mode==1) {
      ## nothing do to
    } else if (mode==2) {
      el[,1:2] <- el[,c(2,1)]
    } else if (mode==3) {
      el <- rbind(el, el[,c(2,1,3)])
    }

    res <- empty_graph(n=n1+n2, directed=directed)
    weight <- list(el[,3])
    names(weight) <- weighted
    res <- add_edges(res, edges=t(as.matrix(el[,1:2])), attr=weight)

  } else {

    if (multiple) {
      el[,3] <- ceiling(el[,3])
      el[,3][ el[,3] < 0 ] <- 0
    } else {
      el[,3] <- el[,3] != 0
    }

    if (!directed || mode==1) {
      ## nothing do to
    } else if (mode==2) {
      el[,1:2] <- el[,c(2,1)]
    } else if (mode==3) {
      el <- rbind(el, el[,c(2,1,3)])
    }

    edges <- unlist(apply(el, 1, function(x) rep(unname(x[1:2]), x[3])))
    res <- graph(n=n1+n2, edges, directed=directed)
  }

  set_vertex_attr(res, "type", value=c(rep(FALSE, n1), rep(TRUE, n2)))
}

graph.incidence.dense <- function(incidence, directed, mode, multiple,
                                  weighted) {

  if (!is.null(weighted)) {
    if (is.logical(weighted) && weighted) {
      weighted <- "weight"
    }
    if (!is.character(weighted)) {
      stop("invalid value supplied for `weighted' argument, please see docs.")
    }

    n1 <- nrow(incidence)
    n2 <- ncol(incidence)
    no.edges <- sum(incidence != 0)
    if (directed && mode==3) { no.edges <- no.edges * 2 }
    edges <- numeric(2*no.edges)
    weight <- numeric(no.edges)
    ptr <- 1
    for (i in seq_len(nrow(incidence))) {
      for (j in seq_len(ncol(incidence))) {
        if (incidence[i,j] != 0) {
          if (!directed || mode==1) {
            edges[2*ptr-1] <- i
            edges[2*ptr] <- n1+j
            weight[ptr] <- incidence[i,j]
            ptr <- ptr + 1
          } else if (mode==2) {
            edges[2*ptr-1] <- n1+j
            edges[2*ptr] <- i
            weight[ptr] <- incidence[i,j]
            ptr <- ptr + 1
          } else if (mode==3) {
            edges[2*ptr-1] <- i
            edges[2*ptr] <- n1+j
            weight[ptr] <- incidence[i,j]
            ptr <- ptr + 1
            edges[2*ptr-1] <- n1+j
            edges[2*ptr] <- i
          }
        }
      }
    }
    res <- empty_graph(n=n1+n2, directed=directed)
    weight <- list(weight)
    names(weight) <- weighted
    res <- add_edges(res, edges, attr=weight)
    res <- set_vertex_attr(res, "type",
                                value=c(rep(FALSE, n1), rep(TRUE, n2)))

  } else {

    mode(incidence) <- "double"
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    ## Function call
    res <- .Call("R_igraph_incidence", incidence, directed, mode, multiple,
                 PACKAGE="igraph")
    res <- set_vertex_attr(res$graph, "type", value=res$types)

  }

  res
}



#' Create graphs from an incidence matrix
#' 
#' \code{graph_from_incidence_matrix} creates a bipartite igraph graph from an incidence
#' matrix.
#' 
#' Bipartite graphs have a \sQuote{\code{type}} vertex attribute in igraph,
#' this is boolean and \code{FALSE} for the vertices of the first kind and
#' \code{TRUE} for vertices of the second kind.
#' 
#' \code{graph_from_incidence_matrix} can operate in two modes, depending on the
#' \code{multiple} argument. If it is \code{FALSE} then a single edge is
#' created for every non-zero element in the incidence matrix. If
#' \code{multiple} is \code{TRUE}, then the matrix elements are rounded up to
#' the closest non-negative integer to get the number of edges to create
#' between a pair of vertices.
#'
#' @aliases graph.incidence
#' @param incidence The input incidence matrix. It can also be a sparse matrix
#' from the \code{Matrix} package.
#' @param directed Logical scalar, whether to create a directed graph.
#' @param mode A character constant, defines the direction of the edges in
#' directed graphs, ignored for undirected graphs. If \sQuote{\code{out}}, then
#' edges go from vertices of the first kind (corresponding to rows in the
#' incidence matrix) to vertices of the second kind (columns in the incidence
#' matrix). If \sQuote{\code{in}}, then the opposite direction is used. If
#' \sQuote{\code{all}} or \sQuote{\code{total}}, then mutual edges are created.
#' @param multiple Logical scalar, specifies how to interpret the matrix
#' elements. See details below.
#' @param weighted This argument specifies whether to create a weighted graph
#' from the incidence matrix. If it is \code{NULL} then an unweighted graph is
#' created and the \code{multiple} argument is used to determine the edges of
#' the graph. If it is a character constant then for every non-zero matrix
#' entry an edge is created and the value of the entry is added as an edge
#' attribute named by the \code{weighted} argument. If it is \code{TRUE} then a
#' weighted graph is created and the name of the edge attribute will be
#' \sQuote{\code{weight}}.
#' @param add.names A character constant, \code{NA} or \code{NULL}.
#' \code{graph_from_incidence_matrix} can add the row and column names of the incidence
#' matrix as vertex attributes. If this argument is \code{NULL} (the default)
#' and the incidence matrix has both row and column names, then these are added
#' as the \sQuote{\code{name}} vertex attribute. If you want a different vertex
#' attribute for this, then give the name of the attributes as a character
#' string. If this argument is \code{NA}, then no vertex attributes (other than
#' type) will be added.
#' @return A bipartite igraph graph. In other words, an igraph graph that has a
#' vertex attribute \code{type}.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{bipartite_graph}} for another way to create bipartite
#' graphs
#' @keywords graphs
#' @examples
#' 
#' inc <- matrix(sample(0:1, 15, repl=TRUE), 3, 5)
#' colnames(inc) <- letters[1:5]
#' rownames(inc) <- LETTERS[1:3]
#' graph_from_incidence_matrix(inc)
#' 
graph_from_incidence_matrix <- function(incidence, directed=FALSE,
                            mode=c("all", "out", "in", "total"),
                            multiple=FALSE, weighted=NULL,
                            add.names=NULL) {
  # Argument checks
  directed <- as.logical(directed)
  mode <- switch(igraph.match.arg(mode), "out"=1, "in"=2, "all"=3, "total"=3)
  multiple <- as.logical(multiple)

  if (inherits(incidence, "Matrix")) {
    res <- graph.incidence.sparse(incidence, directed=directed,
                                  mode=mode, multiple=multiple,
                                  weighted=weighted)
  } else {
    incidence <- as.matrix(incidence)
    res <- graph.incidence.dense(incidence, directed=directed, mode=mode,
                                 multiple=multiple, weighted=weighted)
  }

  ## Add names
  if (is.null(add.names)) {
    if (!is.null(rownames(incidence)) && !is.null(colnames(incidence))) {
      add.names <- "name"
    } else {
      add.names <- NA
    }
  } else if (!is.na(add.names)) {
    if (is.null(rownames(incidence)) || is.null(colnames(incidence))) {
      warning("Cannot add row- and column names, at least one of them is missing")
      add.names <- NA
    }
  }
  if (!is.na(add.names)) {
    res <- set_vertex_attr(res, add.names,
                                value=c(rownames(incidence), colnames(incidence)))
  }
  res
}

#' Create a complete (full) citation graph
#'
#' \code{full_citation_graph} creates a full citation graph. This is a
#' directed graph, where every \code{i->j} edge is present if and only if
#' \eqn{j<i}. If \code{directed=FALSE} then the graph is just a full graph.
#'
#' @aliases graph.full.citation
#' @param n The number of vertices.
#' @param directed Whether to create a directed graph.
#' @return An igraph graph.
#'
#' @family determimistic constructors
#' @export
#' @examples
#' str(full_citation_graph(10))

full_citation_graph <- function(n, directed=TRUE) {
  # Argument checks
  n <- as.integer(n)
  directed <- as.logical(directed)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_full_citation", n, directed,
        PACKAGE="igraph")

  res <- set_graph_attr(res, 'name', 'Full citation graph')
  res
}
