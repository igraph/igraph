
#   IGraph R package
#   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
#   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
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
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
###################################################################

graph <- function( edges, ... ) {
  res <- graph.empty(...)
  if (length(edges) != 0) {

    # Add vertices if needed
    vc <- vcount(res)
    m <- max(edges)
    if (vc < m) {
      res <- add.vertices(res, m-vc)
    }
    
    res <- add.edges(res, edges)
  }
  res
}

graph.star <- function(n, mode="in", center=1, directed=TRUE, ...) {

  if (mode=="out") {
    edges <- c(rep(center, n-1), 2:n)
  } else {
    edges <- c(2:n, rep(center, n-1))
  }

  edges <- as.numeric(t(matrix(edges, nc=2)))
  
  if (mode=="undirected") {
    directed <- FALSE
  }
  
  res <- graph(edges, directed=directed, n=n, ...)
  
  res
}

###################################################################
# Lattices, every kind
###################################################################

graph.lattice <- function(dimvector=NULL,length=NULL, dim=NULL, nei=1,
                          directed=FALSE, mutual=FALSE, circular=FALSE, ...) {

  # Check
  if (is.null(dimvector) && (is.null(length) || is.null(dim))) {
    stop("Either `length' and `dim' or 'dimvector' must be set. See docs.")
  }
  if (!is.null(length) && length < 1) {
    stop("Invalid `length' argument, should be at least one")
  }
  if (!is.null(length) && dim < 1) {
    stop("Invalid `dim' argument, should be at least one")
  }
  if (!is.null(length) && any(dimvector < 1)) {
    stop("Invalid `dimvector', has negative or smaller than one elements")
  }
  if (mutual && !directed) {
    warning("`mutual' specified for undirected graph, proceeding with multiplex edges...")
  }
  if (nei < 1) {
    stop("`nei' should be at least one")
  }
  
  if (!is.null(length)) {
    length <- as.numeric(length)
    dim <- as.numeric(dim)
    dimvector <- rep(length, times=dim)
  } else {
    dimvector <- as.numeric(dimvector)
  }
  nei <- as.numeric(nei)

  n <- prod(dimvector)
  res <- graph.empty(n=n, directed=directed, ...)
  res <- add.edges(res, .Call("REST_create_lattice", dimvector, n,
                              circular, mutual, PACKAGE="igraph"))

  # Connect also to local neighborhood
  if (nei >= 2) {
    neighbors <- lapply(1:length(res), function(a) get.neighborhood(res, a))
    res <- add.edges(res, .Call("REST_connect_neighborhood", neighbors, nei,
                                mutual, PACKAGE="igraph"))
  }
  
  res
}

graph.ring <- function(n, circular=TRUE, directed=FALSE, mutual=FALSE, ...) {
  res <- graph.lattice(dimvector=n, circular=circular, directed=directed,
                       mutual=mutual, ...)
  res
}
