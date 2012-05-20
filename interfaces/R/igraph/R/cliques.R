
#   IGraph R package
#   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_cliques", graph, as.numeric(min), as.numeric(max),
               PACKAGE="igraph")
  lapply(res, function(x) x+1)
}

largest.cliques <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_largest_cliques", graph,
               PACKAGE="igraph")
  lapply(res, function(x) x+1)  
}

maximal.cliques <- function(graph, min=NULL, max=NULL) {
  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }

  if (is.null(min)) { min <- 0 }
  if (is.null(max)) { max <- 0 }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_maximal_cliques", graph, as.numeric(min),
               as.numeric(max),
               PACKAGE="igraph")
  lapply(res, function(x) x+1)
}

clique.number <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
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
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_independent_vertex_sets", graph, as.numeric(min),
               as.numeric(max),
               PACKAGE="igraph")
  lapply(res, function(x) x+1)
}

largest.independent.vertex.sets <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_largest_independent_vertex_sets", graph,
               PACKAGE="igraph")
  lapply(res, function(x) x+1)
}

maximal.independent.vertex.sets <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_maximal_independent_vertex_sets", graph,
               PACKAGE="igraph")
  lapply(res, function(x) x+1)
}

independence.number <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_independence_number", graph,
        PACKAGE="igraph")
}
