
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
#   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
#   02110-1301 USA
#
###################################################################

graph <- function( edges, n=max(edges)+1, directed=TRUE ) {
  .Call("R_igraph_create", as.numeric(edges), as.numeric(n),
        as.logical(directed),
        PACKAGE="igraph")
}

graph.adjacency <- function( adjmatrix, mode="directed" ) {
  mode <- switch(mode,
                 "directed"=0,
                 "undirected"=1,
                 "max"=1,
                 "upper"=2,
                 "lower"=3,
                 "min"=4,
                 "plus"=5)
  attrs <- attributes(adjmatrix)
  adjmatrix <- as.numeric(adjmatrix)
  attributes(adjmatrix) <- attrs
  .Call("R_igraph_graph_adjacency", adjmatrix, as.numeric(mode),
        PACKAGE="igraph")
}
  

graph.star <- function(n, mode="in", center=0 ) {

  if (is.character(mode)) {
    mode <- switch(mode, "out"=0, "in"=1, "undirected"=2)
  }
  .Call("R_igraph_star", as.numeric(n), as.numeric(mode),
        as.numeric(center),
        PACKAGE="igraph")
}

graph.full <- function(n, directed=FALSE, loops=FALSE) {
  .Call("R_igraph_full", as.numeric(n), as.logical(directed),
        as.logical(loops),
        PACKAGE="igraph")
}

###################################################################
# Lattices, every kind
###################################################################

graph.lattice <- function(dimvector=NULL,length=NULL, dim=NULL, nei=1,
                          directed=FALSE, mutual=FALSE, circular=FALSE, ...) {

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
##   res <- graph.empty(n=n, directed=directed, ...)
##   res <- add.edges(res, .Call("REST_create_lattice", dimvector, n,
##                               circular, mutual, PACKAGE="igraph"))

##   # Connect also to local neighborhood
##   if (nei >= 2) {
##     neighbors <- lapply(1:length(res), function(a) get.neighborhood(res, a))
##     res <- add.edges(res, .Call("REST_connect_neighborhood", neighbors, nei,
##                                 mutual, PACKAGE="igraph"))
##   }
  
##   res

  if (is.null(dimvector)) {
    dimvector <- rep(length, dim)
  }
  
  .Call("R_igraph_lattice", as.numeric(dimvector), as.numeric(nei),
        as.logical(directed), as.logical(mutual),
        as.logical(circular),
        PACKAGE="igraph")
}

graph.ring <- function(n, directed=FALSE, mutual=FALSE, circular=TRUE) {
  .Call("R_igraph_ring", as.numeric(n), as.logical(directed),
        as.logical(mutual), as.logical(circular),
        PACKAGE="igraph")
}

###################################################################
# Trees, regular
###################################################################

graph.tree <- function(n, children=2, mode="out") {
  if (is.character(mode)) {
    mode <- switch(mode, "out"=0, "in"=1, "undirected"=2);
  }

  .Call("R_igraph_tree", as.numeric(n), as.numeric(children),
        as.numeric(mode),
        PACKAGE="igraph")
}

###################################################################
# The graph atlas
###################################################################

graph.atlas <- function(n) {

  .Call("R_igraph_atlas", as.numeric(n),
        PACKAGE="igraph")
}
