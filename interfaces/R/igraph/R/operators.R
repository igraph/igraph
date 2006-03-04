
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

graph.disjoint.union <- function(...) {

  graphs <- unlist(recursive=FALSE, lapply(list(...), function(l) {
    if (is.igraph(l)) list(l) else l
  } ))
  
  .Call("R_igraph_disjoint_union", graphs,
        PACKAGE="igraph")
}

"%du%" <- function(x,y) {
  graph.disjoint.union(x,y)
}

graph.union <- function(...) {

  graphs <- unlist(recursive=FALSE, lapply(list(...), function(l) {
    if (is.igraph(l)) list(l) else l
  } ))
  
  .Call("R_igraph_union", graphs,
        PACKAGE="igraph")
}

"%u%" <- function(x,y) {
  graph.union(x,y)
}

graph.intersection <- function(...) {

  graphs <- unlist(recursive=FALSE, lapply(list(...), function(l) {
    if (is.igraph(l)) list(l) else l
  } ))
  
  .Call("R_igraph_intersection", graphs,
        PACKAGE="igraph")
}

"%s%" <- function(x,y) {
  graph.intersection(x,y)
}

graph.difference <- function(big, small) {

  if (!is.igraph(big) || !is.igraph(small)) {
    stop("argument is not a graph")
  }
  
  .Call("R_igraph_difference", big, small,
        PACKAGE="igraph")
}
    
"%m%" <- function(x,y) {
  graph.difference(x,y)
}

graph.complementer <- function(graph, loops=FALSE) {

  .Call("R_igraph_complementer", graph, as.logical(loops),
        PACKAGE="igraph")
}

graph.compose <- function(g1, g2) {

  .Call("R_igraph_compose", g1, g2,
        PACKAGE="igraph")
}

"%c%" <- function(x,y) {
  graph.compose(x,y)
}
