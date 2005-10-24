
#   IGraph R package
#   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
###################################################################

ii.create <- function(graph, type, vid=NULL, mode=NULL) {
  if (type=="vid") {
    .Call("R_igraph_iterator_vid", graph, PACKAGE="igraph")
  } else if (type=="eid") {
    .Call("R_igraph_iterator_eid", graph, PACKAGE="igraph")
  } else if (type=="efromorder") {
    .Call("R_igraph_iterator_efromorder", graph, PACKAGE="igraph")
  } else if (type=="eneis") {
    if (is.character(mode)) {
      mode <- switch(mode, "out"=1, "in"=2, "all"=3, "total"=3)
    }
    .Call("R_igraph_iterator_eneis", graph, as.numeric(vid),
          as.numeric(mode),
          PACKAGE="igraph")
  } else if (type=="vneis") {
    if (is.character(mode)) {
      mode <- switch(mode, "out"=1, "in"=2, "all"=3, "total"=3)
    }
    .Call("R_igraph_iterator_vneis", graph, as.numeric(vid),
          as.numeric(mode),
          PACKAGE="igraph")
  }
}

ii.next <- function(graph, iterator) {
  .Call("R_igraph_iterator_next", graph, iterator, PACKAGE="igraph")
}

ii.prev <- function(graph, iterator) {
  .Call("R_igraph_iterator_prev", graph, iterator, PACKAGE="igraph")
}

ii.end <- function(graph, iterator) {
  .Call("R_igraph_iterator_end", graph, iterator, PACKAGE="igraph")
}

ii.get.vertex.nei <- function(graph, iterator) {
  .Call("R_igraph_iterator_get_vertex_nei", graph, iterator, PACKAGE="igraph")
}

ii.reset <- function(graph, iterator) {
  .Call("R_igraph_iterator_reset", graph, iterator, PACKAGE="igraph")
}

ii.get.vertex <- function(graph, iterator) {
  .Call("R_igraph_iterator_get_vertex", graph, iterator, PACKAGE="igraph")
}

ii.get.from <- function(graph, iterator) {
  .Call("R_igraph_iterator_from", graph, iterator, PACKAGE="igraph")
}

ii.get.to <- function(graph, iterator) {
  .Call("R_igraph_iterator_get_to", graph, iterator, PACKAGE="igraph")
}

ii.get.edge <- function(graph, iterator) {
  .Call("R_igraph_iterator_edge", graph, iterator, PACKAGE="igraph")
}
