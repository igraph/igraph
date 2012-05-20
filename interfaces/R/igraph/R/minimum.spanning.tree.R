
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

minimum.spanning.tree <- function(graph, weights=NULL,
                                  algorithm=NULL, ...) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(algorithm)) {
    if (!is.null(weights) || "weight" %in% list.edge.attributes(graph)) {
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
    if (is.null(weights) && ! "weight" %in% list.edge.attributes(graph)) {
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
