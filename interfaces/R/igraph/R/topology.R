
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

graph.isoclass <- function(graph) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_isoclass_34", graph,
        PACKAGE="igraph")
}

graph.isomorphic <- function(graph1, graph2) {

  if (!is.igraph(graph1) || !is.igraph(graph2)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_isomorphic_34", graph1, graph2,
        PACKAGE="igraph")
}

graph.isocreate <- function(size, number, directed=TRUE) {

  .Call("R_igraph_isoclass_create", as.numeric(size),
        as.numeric(number), as.logical(directed),
        PACKAGE="igraph")
}

graph.isomorphic.vf2 <- function(graph1, graph2) {

  if (!is.igraph(graph1) || !is.igraph(graph2)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_isomorphic_vf2", graph1, graph2,
        PACKAGE="igraph")
}  

