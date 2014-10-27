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

#' @export

arpack_defaults <- list(bmat="I", n=0, which="XX", nev=1, tol=0.0,
                              ncv=3, ldv=0, ishift=1, maxiter=3000, nb=1,
                              mode=1, start=0, sigma=0.0, sigmai=0.0)

#' @export

arpack <- function(func, extra=NULL, sym=FALSE, options=arpack_defaults,
                   env=parent.frame(), complex=!sym) {

  if (!is.list(options) ||
      (is.null(names(options)) && length(options) != 0)) {
    stop("options must be a named list")
  }
  if (any(names(options) == "")) {
    stop("all options must be named")
  }
  if (any(! names(options) %in% names(arpack_defaults))) {
    stop("unkown ARPACK option(s): ",
         paste(setdiff(names(options), names(arpack_defaults)),
                       collapse=", "))
  }
  
  options.tmp <- arpack_defaults
  options.tmp[ names(options) ] <- options
  options <- options.tmp

  if (sym && complex) {
    complex <- FALSE
    warning("Symmetric matrix, setting `complex' to FALSE")
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_arpack", func, extra, options, env, sym,
               PACKAGE="igraph")

  if (complex) {
    rew <- arpack.unpack.complex(res$vectors, res$values,
                                 min(res$options$nev, res$options$nconv))
    res$vectors <- rew$vectors
    res$values <- rew$values

    res$values <- apply(res$values, 1, function(x) x[1]+x[2]*1i)
    dim(res$vectors) <- c(nrow(res$vectors)*2, ncol(res$vectors)/2)
    res$vectors <- apply(res$vectors, 2, function(x) {
      l <- length(x)/2
      x[1:l] + x[(l+1):length(x)]*1i
    })
  } else {
    if (is.matrix(res$values)) {
      if (!all(res$values[,2]==0)) {
        warning("Dropping imaginary parts of eigenvalues")
      }
      res$values <- res$values[,1]
    }
    res$vectors <- res$vectors[,1:length(res$values)]
  }
  
  res
}

arpack.unpack.complex <- function(vectors, values, nev) {
  # Argument checks
  vectors <- as.matrix(structure(as.double(vectors), dim=dim(vectors)))
  values <- as.matrix(structure(as.double(values), dim=dim(values)))
  nev <- as.integer(nev)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_arpack_unpack_complex", vectors, values, nev,
        PACKAGE="igraph")

  res
}



#' Find subgraph centrality scores of network positions
#' 
#' Subgraph centrality of a vertex measures the number of subgraphs a vertex
#' participates in, weighting them according to their size.
#' 
#' The subgraph centrality of a vertex is defined as the number of closed loops
#' originating at the vertex, where longer loops are exponentially
#' downweighted.
#' 
#' Currently the calculation is performed by explicitly calculating all
#' eigenvalues and eigenvectors of the adjacency matrix of the graph. This
#' effectively means that the measure can only be calculated for small graphs.
#'
#' @aliases subgraph.centrality
#' @param graph The input graph, it should be undirected, but the
#' implementation does not check this currently.
#' @param diag Boolean scalar, whether to include the diagonal of the adjacency
#' matrix in the analysis. Giving \code{FALSE} here effectively eliminates the
#' loops edges from the graph before the calculation.
#' @return A numeric vector, the subgraph centrality scores of the vertices.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com} based on the Matlab
#' code by Ernesto Estrada
#' @seealso \code{\link{eigen_centrality}}, \code{\link{page_rank}}
#' @references Ernesto Estrada, Juan A. Rodriguez-Velazquez: Subgraph
#' centrality in Complex Networks. \emph{Physical Review E} 71, 056103 (2005).
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- sample_pa(100, m=4, dir=FALSE)
#' sc <- subgraph_centrality(g)
#' cor(degree(g), sc)
#' 
subgraph_centrality <- function(graph, diag=FALSE) {
  A <- as_adj(graph)
  if (!diag) { diag(A) <- 0 }
  eig <- eigen(A)
  res <- as.vector(eig$vectors^2 %*% exp(eig$values))
  if (igraph_opt("add.vertex.names") && is_named(graph)) { 
    names(res) <- vertex_attr(graph, "name") 
  }
  res
}

#' @export

eigen_defaults <- list(pos="LM", howmany=1L, il=-1L, iu=-1L,
                             vl=-Inf, vu=Inf, vestimate=0L,
                             balance="none")
