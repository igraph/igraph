
#   IGraph R package
#   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

bipartite.projection <- function(graph, types=NULL,
                                 multiplicity=TRUE, probe1=NULL) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  if (is.null(types) && "type" %in% list.vertex.attributes(graph)) { 
  types <- V(graph)$type 
  } 
  if (!is.null(types)) { 
  types <- as.logical(types) 
  } else { 
  stop("Not a bipartite graph, supply `types' argument") 
  }
  if (!is.null(probe1)) {
    probe1 <- as.igraph.vs(graph, probe1)-1
  } else {
    probe1 <- -1
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_bipartite_projection", graph, types, probe1,
        PACKAGE="igraph")
  if (multiplicity) {
    E(res[[1]])$weight <- res[[3]]
    E(res[[2]])$weight <- res[[4]]
  }
  res[1:2]
}
