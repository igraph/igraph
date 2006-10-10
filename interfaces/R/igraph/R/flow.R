
#   IGraph R package
#   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
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

graph.maxflow <- function(graph, source, target, capacity=NULL) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (is.null(capacity)) {
    if ("capacity" %in% list.edge.attributes(graph)) {
      capacity <- E(graph)$capacity
    } else {
      stop("capacity argument is not given and no `capacity' attribute")
    }
  }
  capacity <- as.numeric(capacity)
  
  .Call("R_igraph_maxflow", graph, as.numeric(source), as.numeric(target),
        capacity,
        PACKAGE="igraph")
}

graph.mincut <- function(graph, source=NULL, target=NULL, capacity=NULL) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (is.null(capacity)) {
    if ("capacity" %in% list.edge.attributes(graph)) {
      capacity <- E(graph)$capacity
    } else {
      stop("capacity argument is not given and no `capacity' attribute")
    }
  }
  if (is.null(source) && !is.null(target) ||
      is.null(target) && !is.null(source)) {
    stop("Please give both source and target oe neither")
  }
  capacity <- as.numeric(capacity)
  
  if (is.null(target) && is.null(source)) {
    .Call("R_igraph_mincut_value", graph, capacity,
          PACKAGE="igraph")
  } else {
    .Call("R_igraph_st_mincut_value", graph, as.numeric(source),
          as.numeric(target), capacity,
          PACKAGE="igraph")
  }
}

vertex.connectivity <- function(graph, source=NULL, target=NULL) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(source) && is.null(target)) {
    .Call("R_igraph_vertex_connectivity", graph,
          PACKAGE="igraph")
  } else if (!is.null(source) && !is.null(target)) {
    .Call("R_igraph_st_vertex_connectivity", graph, as.numeric(source),
          as.numeric(target),
          PACKAGE="igraph")
  } else {
    stop("either give both source and target or neither")
  }
}

edge.connectivity <- function(graph, source=NULL, target=NULL) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(source) && is.null(target)) {    
    .Call("R_igraph_edge_connectivity", graph,
          PACKAGE="igraph")
  } else if (!is.null(source) && !is.null(target)) {
    .Call("R_igraph_st_edge_connectivity", graph,
          as.numeric(source), as.numeric(target),
          PACKAGE="igraph")
  } else {
    stop("either give both source and target or neither")
  }
}

edge.disjoint.paths <- function(graph, source, target) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  .Call("R_igraph_edge_disjoint_paths", graph,
        as.numeric(source), as.numeric(target),
        PACKAGE="igraph")
}

vertex.disjoint.paths <- function(graph, source=NULL, target=NULL) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  .Call("R_igraph_vertex_disjoint_paths", graph, as.numeric(source),
        as.numeric(target),
        PACKAGE="igraph")
}

graph.adhesion <- function(graph) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  
  .Call("R_igraph_adhesion", graph,
        PACKAGE="igraph")
}

graph.cohesion <- function(graph) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  .Call("R_igraph_cohesion", graph,
        PACKAGE="igraph")
}

