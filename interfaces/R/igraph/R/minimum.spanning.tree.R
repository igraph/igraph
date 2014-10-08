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



#' Minimum spanning tree
#' 
#' A subgraph of a connected graph is a \emph{minimum spanning tree} if it is
#' tree, and the sum of its edge weights are the minimal among all tree
#' subgraphs of the graph. A minimum spanning forest of a graph is the graph
#' consisting of the minimum spanning trees of its components.
#' 
#' If the graph is unconnected a minimum spanning forest is returned.
#'
#' @aliases minimum.spanning.tree
#' @param graph The graph object to analyze.
#' @param weights Numeric algorithm giving the weights of the edges in the
#' graph. The order is determined by the edge ids. This is ignored if the
#' \code{unweighted} algorithm is chosen
#' @param algorithm The algorithm to use for calculation. \code{unweighted} can
#' be used for unwieghted graphs, and \code{prim} runs Prim's algorithm for
#' weighted graphs.  If this is \code{NULL} then igraph tries to select the
#' algorithm automatically: if the graph has an edge attribute called
#' \code{weight} of the \code{weights} argument is not \code{NULL} then Prim's
#' algorithm is chosen, otherwise the unwweighted algorithm is performed.
#' @param \dots Additional arguments, unused.
#' @return A graph object with the minimum spanning forest. (To check that it
#' is a tree check that the number of its edges is \code{vcount(graph)-1}.)
#' The edge and vertex attributes of the original graph are preserved in the
#' result.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{components}}
#' @references Prim, R.C. 1957. Shortest connection networks and some
#' generalizations \emph{Bell System Technical Journal}, 37 1389--1401.
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- sample_gnp(100, 3/100)
#' g_mst <- mst(g)
#' 
mst <- function(graph, weights=NULL,
                                  algorithm=NULL, ...) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(algorithm)) {
    if (!is.null(weights) || "weight" %in% edge_attr_names(graph)) {
      algorithm <- "prim"
    } else {
      algorithm <- "unweighted"
    }
  }
  
  if (algorithm=="unweighted") {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    .Call("R_igraph_minimum_spanning_tree_unweighted", graph,
          PACKAGE="igraph")
  } else if (algorithm=="prim") {
    if (is.null(weights) && ! "weight" %in% edge_attr_names(graph)) {
      stop("edges weights must be supplied for Prim's algorithm")
    } else if (is.null(weights)) {
      weights <- E(graph)$weight
    }
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    .Call("R_igraph_minimum_spanning_tree_prim", graph, as.numeric(weights),
          PACKAGE="igraph")    
  } else {
    stop("Invalid algorithm")
  }
}
