
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

  oclass <- class(graph)
  graph <- unclass(graph)
  for (i in seq(attrs)) {
    graph[[9]][[4]][[nam[i]]][idx] <- attrs[[nam[i]]]
  }  
  class(graph) <- oclass
  
  graph
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

  oclass <- class(graph)
  graph <- unclass(graph)
  for (i in seq(attrs)) {
    graph[[9]][[3]][[nam[i]]][idx] <- attrs[[nam[i]]]
  }
  class(graph) <- oclass
                  
  graph
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
  if (is.character(mode)) {
    mode <- switch(mode, "out"=1, "in"=2, "all"=3, "total"=3)
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
  matrix(res, nc=2, byrow=TRUE)+1
}

get.edge.ids <- function(graph, vp, directed=TRUE, error=FALSE, multi=FALSE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_get_eids", graph, as.igraph.vs(graph, vp)-1,
        as.logical(directed), as.logical(error), as.logical(multi),
        PACKAGE="igraph")+1
}  
