
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

cliques <- function(graph, min=NULL, max=NULL) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(min)) {
    min <- 0
  }
  if (is.null(max)) {
    max <- 0
  }

  .Call("R_igraph_cliques", graph, as.numeric(min), as.numeric(max),
        PACKAGE="igraph")
}

largest.cliques <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  .Call("R_igraph_largest_cliques", graph,
        PACKAGE="igraph")
}

maximal.cliques <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }

  .Call("R_igraph_maximal_cliques", graph,
        PACKAGE="igraph")
}

clique.number <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }
  
  .Call("R_igraph_clique_number", graph,
        PACKAGE="igraph")
}

independent.vertex.sets <- function(graph, min=NULL, max=NULL) {
  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }

  if (is.null(min)) {
    min <- 0
  }

  if (is.null(max)) {
    max <- 0
  }
  
  .Call("R_igraph_independent_vertex_sets", graph, as.numeric(min),
        as.numeric(max),
        PACKAGE="igraph")
}

largest.independent.vertex.sets <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }

  .Call("R_igraph_largest_independent_vertex_sets", graph,
        PACKAGE="igraph")
}

maximal.independent.vertex.sets <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }

  .Call("R_igraph_maximal_independent_vertex_sets", graph,
        PACKAGE="igraph")
}

independence.number <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }

  .Call("R_igraph_independence_number", graph,
        PACKAGE="igraph")
}
