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



#' Cocitation coupling
#' 
#' Two vertices are cocited if there is another vertex citing both of them.
#' \code{cocitation} siply counts how many types two vertices are cocited. The
#' bibliographic coupling of two vertices is the number of other vertices they
#' both cite, \code{bibcoupling} calculates this.
#' 
#' \code{cocitation} calculates the cocitation counts for the vertices in the
#' \code{v} argument and all vertices in the graph.
#' 
#' \code{bibcoupling} calculates the bibliographic coupling for vertices in
#' \code{v} and all vertices in the graph.
#' 
#' Calculating the cocitation or bibliographic coupling for only one vertex
#' costs the same amount of computation as for all vertices. This might change
#' in the future.
#' 
#' @aliases cocitation bibcoupling
#' @param graph The graph object to analyze
#' @param v Vertex sequence or numeric vector, the vertex ids for which the
#' cocitation or bibliographic coupling values we want to calculate. The
#' default is all vertices.
#' @return A numeric matrix with \code{length(v)} lines and
#' \code{vcount(graph)} columns. Element \code{(i,j)} contains the cocitation
#' or bibliographic coupling for vertices \code{v[i]} and \code{j}.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- make_ring(10)
#' cocitation(g)
#' bibcoupling(g)
#' 
cocitation <- function(graph, v=V(graph)) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  v <- as.igraph.vs(graph, v)
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_cocitation", graph, v-1,
               PACKAGE="igraph")
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    rownames(res) <- vertex_attr(graph, "name", v)
    colnames(res) <- vertex_attr(graph, "name")
  }
  res
}

#' @export

bibcoupling <- function(graph, v=V(graph)) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  v <- as.igraph.vs(graph, v)
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_bibcoupling", graph, v-1,
               PACKAGE="igraph")
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    rownames(res) <- vertex_attr(graph, "name", v)
    colnames(res) <- vertex_attr(graph, "name")
  }
  res
}
