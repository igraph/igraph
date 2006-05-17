
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

graph.empty <- function(n=0, directed=TRUE) {
  .Call("R_igraph_empty", as.numeric(n), as.logical(directed),
        PACKAGE="igraph")
}

add.edges <- function(graph, edges) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_add_edges", graph, as.numeric(edges),
        PACKAGE="igraph")
}

add.vertices <- function(graph, nv) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_add_vertices", graph, as.numeric(nv),
        PACKAGE="igraph")
}

delete.edges <- function(graph, edges) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_delete_edges", graph, as.numeric(edges),
        PACKAGE="igraph")
}

delete.vertices <- function(graph, v) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_delete_vertices", graph, as.igraph.vs(v),
        PACKAGE="igraph")
}

###################################################################
# Structure query
###################################################################
  
vcount <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_vcount", graph,
        PACKAGE="igraph")
}
  
ecount <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
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
  .Call("R_igraph_neighbors", graph, as.numeric(v),
        as.numeric(mode),
        PACKAGE="igraph")
}

is.directed <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_is_directed", graph,
        PACKAGE="igraph")
}
