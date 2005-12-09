
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

igraph.vs.all <- function(graph) {
  .Call("R_igraph_vs_all", graph, PACKAGE="igraph")
}

igraph.es.all <- function(graph) {
  .Call("R_igraph_es_all", graph, PACKAGE="igraph")
}

igraph.es.fromorder <- function(graph) {
  .Call("R_igraph_es_fromorder", graph, PACKAGE="igraph")
}

igraph.es.adj <- function(graph, vid, mode="all") {
  if (is.character(mode)) {
    mode <- switch(mode, "out"=1, "in"=2, "all"=3, "total"=3)
  }
  .Call("R_igraph_es_adj", graph, as.numeric(vid), as.numeric(mode),
        PACKAGE="igraph")
}

igraph.vs.adj <- function(graph, vid, mode="all") {
  if (is.character(mode)) {
    mode <- switch(mode, "out"=1, "in"=2, "all"=3, "total"=3)
  }
  .Call("R_igraph_vs_adj", graph, as.numeric(vid), as.numeric(mode),
        PACKAGE="igraph")
}

igraph.vs.vector <- function(graph, v) {
  .Call("R_igraph_vs_vector", graph, as.numeric(v),
        PACKAGE="igraph")
}

igraph.es.vector <- function(graph, v) {
  .Call("R_igraph_es_vector", graph, as.numeric(v),
        PACKAGE="igraph")
}

igraph.vs.next <- function(graph, iterator) {
  .Call("R_igraph_vs_next", graph, iterator, PACKAGE="igraph")
}

igraph.vs.end <- function(graph, iterator) {
  .Call("R_igraph_vs_end", graph, iterator, PACKAGE="igraph")
}

igraph.vs.get <- function(graph, iterator) {
  .Call("R_igraph_vs_get", graph, iterator, PACKAGE="igraph")
}

igraph.vs.reset <- function(graph, iterator) {
  .Call("R_igraph_vs_reset", graph, iterator, PACKAGE="igraph")
}

igraph.es.next <- function(graph, iterator) {
  .Call("R_igraph_es_next", graph, iterator, PACKAGE="igraph")
}

igraph.es.end <- function(graph, iterator) {
  .Call("R_igraph_es_end", graph, iterator, PACKAGE="igraph")
}

igraph.es.get <- function(graph, iterator) {
  .Call("R_igraph_es_get", graph, iterator, PACKAGE="igraph")
}

igraph.es.reset <- function(graph, iterator) {
  .Call("R_igraph_es_reset", graph, iterator, PACKAGE="igraph")
}

igraph.es.from <- function(graph, iterator) {
  .Call("R_igraph_es_from", graph, iterator, PACKAGE="igraph")
}

igraph.es.to <- function(graph, iterator) {
  .Call("R_igraph_es_to", graph, iterator, PACKAGE="igraph")
}

#####
# internal helper functions

as.igraph.vs <- function(g, it) {
  if (class(it) == "igraph.vs") {
    it
  } else if (is.numeric(it)) {
    igraph.vs.vector(g, it)
  } else {
    stop("Cannot interpret vertex set")
  }
}
