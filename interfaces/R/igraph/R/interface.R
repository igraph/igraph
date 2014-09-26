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

###################################################################
# Structure building
###################################################################

add.edges <- function(graph, edges, ..., attr=list()) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  attrs <- list(...)
  attrs <- append(attrs, attr)
  nam <- names(attrs)
  if (length(attrs) != 0 && (is.null(nam) || any(nam==""))) {
    stop("please supply names for attributes")
  }
  
  edges.orig <- ecount(graph)  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  graph <- .Call("R_igraph_add_edges", graph, as.igraph.vs(graph, edges)-1,
                 PACKAGE="igraph")
  edges.new <- ecount(graph)
  
  if (edges.new-edges.orig != 0) {
    idx <- seq(edges.orig+1, edges.new)
  } else {
    idx <- numeric()
  }

  eattrs <- .Call("R_igraph_mybracket2", graph, 9L, 4L, PACKAGE="igraph")
  for (i in seq(attrs)) { eattrs[[nam[i]]][idx] <- attrs[[nam[i]]] }

  .Call("R_igraph_mybracket2_set", graph, 9L, 4L, eattrs, PACKAGE="igraph")
}

add.vertices <- function(graph, nv, ..., attr=list()) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  attrs <- list(...)
  attrs <- append(attrs, attr)
  nam <- names(attrs)
  if (length(attrs) != 0 && (is.null(nam) || any(nam==""))) {
    stop("please supply names for attributes")
  }

  vertices.orig <- vcount(graph)  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  graph <- .Call("R_igraph_add_vertices", graph, as.numeric(nv),
                 PACKAGE="igraph")
  vertices.new <- vcount(graph)

  if (vertices.new-vertices.orig != 0) {
    idx <- seq(vertices.orig+1, vertices.new)
  } else {
    idx <- numeric()
  }

  vattrs <- .Call("R_igraph_mybracket2", graph, 9L, 3L, PACKAGE="igraph")
  for (i in seq(attrs)) { vattrs[[nam[i]]][idx] <- attrs[[nam[i]]] }

  .Call("R_igraph_mybracket2_set", graph, 9L, 3L, vattrs, PACKAGE="igraph")
}

delete.edges <- function(graph, edges) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_delete_edges", graph, as.igraph.es(graph, edges)-1,
        PACKAGE="igraph")
}

delete.vertices <- function(graph, v) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_delete_vertices", graph, as.igraph.vs(graph, v)-1,
        PACKAGE="igraph")
}

###################################################################
# Structure query
###################################################################
  
ecount <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_ecount", graph,
        PACKAGE="igraph")
}
 
neighbors <- function(graph, v, mode=1) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (is.character(mode)) {
    mode <- switch(mode, "out"=1, "in"=2, "all"=3, "total"=3)
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_neighbors", graph, as.igraph.vs(graph, v)-1,
               as.numeric(mode),
               PACKAGE="igraph")
  res+1
}

incident <- function(graph, v, mode=c("all", "out", "in", "total")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (is.directed(graph)) {
    mode <- igraph.match.arg(mode)
    mode <- switch(mode, "out"=1, "in"=2, "all"=3, "total"=3)
  } else {
    mode=1
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_incident", graph, as.igraph.vs(graph, v)-1,
               as.numeric(mode),
               PACKAGE="igraph")
  res+1
}  

is.directed <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_is_directed", graph,
        PACKAGE="igraph")
}

get.edges <- function(graph, es) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_edges", graph, as.igraph.es(graph, es)-1,
               PACKAGE="igraph")
  matrix(res, ncol=2, byrow=TRUE)+1
}



#' Find the edge ids based on the incident vertices of the edges
#' 
#' Find the edges in an igraph graph that have the specified end points. This
#' function handles multi-graph (graphs with multiple edges) and can consider
#' or ignore the edge directions in directed graphs.
#' 
#' igraph vertex ids are natural numbers, starting from one, up to the number
#' of vertices in the graph. Similarly, edges are also numbered from one, up to
#' the number of edges.
#' 
#' This function allows finding the edges of the graph, via their incident
#' vertices.
#' 
#' @param graph The input graph.
#' @param vp The indicent vertices, given as vertex ids or symbolic vertex
#' names. They are interpreted pairwise, i.e. the first and second are used for
#' the first edge, the third and fourth for the second, etc.
#' @param directed Logical scalar, whether to consider edge directions in
#' directed graphs. This argument is ignored for undirected graphs.
#' @param error Logical scalar, whether to report an error if an edge is not
#' found in the graph. If \code{FALSE}, then no error is reported, and zero is
#' returned for the non-existant edge(s).
#' @param multi Logical scalar, whether to handle multiple edges properly. If
#' \code{FALSE}, and a pair of vertices are given twice (or more), then always
#' the same edge id is reported back for them. If \code{TRUE}, then the edge
#' ids of multiple edges are correctly reported.
#' @return A numeric vector of edge ids, one for each pair of input vertices.
#' If there is no edge in the input graph for a given pair of vertices, then
#' zero is reported. (If the \code{error} argument is \code{FALSE}.)
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @keywords graphs
#' @examples
#' 
#' g <- graph.ring(10)
#' ei <- get.edge.ids(g, c(1,2, 4,5))
#' E(g)[ei]
#' 
#' ## non-existant edge
#' get.edge.ids(g, c(2,1, 1,4, 5,4))
#' 
get.edge.ids <- function(graph, vp, directed=TRUE, error=FALSE, multi=FALSE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_get_eids", graph, as.igraph.vs(graph, vp)-1,
        as.logical(directed), as.logical(error), as.logical(multi),
        PACKAGE="igraph")+1
}  
