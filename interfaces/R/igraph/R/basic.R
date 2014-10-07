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



#' Is this object a graph?
#' 
#' \code{is.graph} makes its decision based on the class attribute of the
#' object.
#' 
#' @aliases is.igraph
#' @param graph An R object.
#' @return A logical constant, \code{TRUE} if argument \code{graph} is a graph
#' object.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @keywords graphs
#' @examples
#' 
#' g <- ring(10)
#' is_igraph(g)
#' is_igraph(numeric(10))
#' 
is_igraph <- function(graph){
  "igraph" %in% class(graph)
}

is_directed <- function(graph) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_is_directed", graph,
        PACKAGE="igraph")
}

get.edge <- function(graph, id) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  id <- as.numeric(id)
  ec <- ecount(graph)
  
  if (id < 1 || id > ec) {
    stop("No such edge")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_get_edge", graph, as.numeric(id)-1,
               PACKAGE="igraph")
  res+1
}


