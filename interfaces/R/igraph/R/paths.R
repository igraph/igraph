## -----------------------------------------------------------------------
##
##   IGraph R package
##   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
##   334 Harvard street, Cambridge, MA 02139 USA
##
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
##   02110-1301 USA
##
## -----------------------------------------------------------------------

#' List all simple paths from one source
#' 
#' This function lists are simple paths from one source vertex to another
#' vertex or vertices. A path is simple if the vertices it visits are not
#' visited more than once.
#' 
#' Note that potentially there are exponentially many paths between two
#' vertices of a graph, and you may run out of memory when using this
#' function, if your graph is lattice-like.
#' 
#' This function currently ignored multiple and loop edges.
#' 
#' @param graph The input graph.
#' @param from The source vertex.
#' @param to The target vertex of vertices. Defaults to all vertices.
#' @param mode Character constant, gives whether the shortest paths to or
#'   from the given vertices should be calculated for directed graphs. If
#'   \code{out} then the shortest paths \emph{from} the vertex, if \code{in}
#'   then \emph{to} it will be considered. If \code{all}, the default, then
#'   the corresponding undirected graph will be used, ie. not directed paths
#'   are searched. This argument is ignored for undirected graphs.
#' @return A list of integer vectors, each integer vector is a path from
#'   the source vertex to one of the target vertices. A path is given by its
#'   vertex ids.
#' @keywords graphs
#' @examples
#' 
#' g <- make_ring(10)
#' all_simple_paths(g, 1, 5)
#' all_simple_paths(g, 1, c(3,5))
#' 
#' @export

all_simple_paths <- function(graph, from, to = V(graph),
                             mode = c("out", "in", "all", "total")) {
  ## Argument checks
  if (!is_igraph(graph)) stop("Not a graph object")
  from <- as.igraph.vs(graph, from)
  to <- as.igraph.vs(graph, to)
  mode <- switch(igraph.match.arg(mode), "out" = 1, "in" = 2, "all" = 3,
                 "total" = 3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )

  ## Function call
  res <- .Call("R_igraph_get_all_simple_paths", graph, from - 1, to - 1,
                mode, PACKAGE = "igraph")
  res <- get.all.simple.paths.pp(res)

  if (igraph_opt("return.vs.es")) { 
    res <- lapply(res, create_vs, graph = graph)
  }
  res
}
