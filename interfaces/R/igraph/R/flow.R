
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
    }
  }
  if (!is.null(capacity)) {
    capacity <- as.numeric(capacity)
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_maxflow", graph, as.igraph.vs(graph, source),
        as.igraph.vs(graph, target), capacity,
        PACKAGE="igraph0")
}

graph.mincut <- function(graph, source=NULL, target=NULL, capacity=NULL,
                         value.only=TRUE) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (is.null(capacity)) {
    if ("capacity" %in% list.edge.attributes(graph)) {
      capacity <- E(graph)$capacity
    }
  }
  if (is.null(source) && !is.null(target) ||
      is.null(target) && !is.null(source)) {
    stop("Please give both source and target or neither")
  }
  if (!is.null(capacity)) {
    capacity <- as.numeric(capacity)
  }

  value.only <- as.logical(value.only)
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  if (is.null(target) && is.null(source)) {
    if (value.only) {
      .Call("R_igraph_mincut_value", graph, capacity,
            PACKAGE="igraph0")
    } else {
      .Call("R_igraph_mincut", graph, capacity,
            PACKAGE="igraph0")
    }
  } else {
    if (value.only) {
      .Call("R_igraph_st_mincut_value", graph, as.igraph.vs(graph, source),
            as.igraph.vs(graph, target), capacity,
            PACKAGE="igraph0")
    } else {
      stop("Calculating minimum s-t cuts is not implemented yet")
    }
  }
}

vertex.connectivity <- function(graph, source=NULL, target=NULL, checks=TRUE) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(source) && is.null(target)) {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
    .Call("R_igraph_vertex_connectivity", graph, as.logical(checks),
          PACKAGE="igraph0")
  } else if (!is.null(source) && !is.null(target)) {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
    .Call("R_igraph_st_vertex_connectivity", graph, as.igraph.vs(graph, source),
          as.igraph.vs(graph, target),
          PACKAGE="igraph0")
  } else {
    stop("either give both source and target or neither")
  }
}

edge.connectivity <- function(graph, source=NULL, target=NULL, checks=TRUE) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(source) && is.null(target)) {    
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
    .Call("R_igraph_edge_connectivity", graph, as.logical(checks),
          PACKAGE="igraph0")
  } else if (!is.null(source) && !is.null(target)) {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
    .Call("R_igraph_st_edge_connectivity", graph,
          as.igraph.vs(graph, source), as.igraph.vs(graph, target),
          PACKAGE="igraph0")
  } else {
    stop("either give both source and target or neither")
  }
}

edge.disjoint.paths <- function(graph, source, target) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_edge_disjoint_paths", graph,
        as.igraph.vs(graph, source), as.igraph.vs(graph, target),
        PACKAGE="igraph0")
}

vertex.disjoint.paths <- function(graph, source=NULL, target=NULL) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_vertex_disjoint_paths", graph, as.igraph.vs(graph, source),
        as.igraph.vs(graph, target),
        PACKAGE="igraph0")
}

graph.adhesion <- function(graph, checks=TRUE) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_adhesion", graph, as.logical(checks),
        PACKAGE="igraph0")
}

graph.cohesion <- function(graph, checks=TRUE) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_cohesion", graph, as.logical(checks),
        PACKAGE="igraph0")
}

