
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

###################################################################
# Structural properties
###################################################################

diameter <- function(graph, directed=TRUE, unconnected=TRUE) {
  .Call("R_igraph_diameter", graph, as.logical(directed),
        as.logical(unconnected),
        PACKAGE="igraph")
}

average.path.length <- function(graph, directed=TRUE, unconnected=TRUE) {
  .Call("R_igraph_average_path_length", graph, as.logical(directed),
        as.logical(unconnected),
        PACKAGE="igraph")
}

degree <- function(graph, v=igraph.vs.all(graph),
                   mode="total", loops=TRUE){
  if (is.character(mode)) {
    mode <- switch(mode, "out"=1, "in"=2, "all"=3, "total"=3)
  }
  
  .Call("R_igraph_degree", graph, as.igraph.vs(graph, v), as.numeric(mode),
        as.logical(loops), PACKAGE="igraph")
}
  
degree.distribution <- function(graph, cumulative=FALSE, ...) {
  
  cs <- degree(graph, ...)
  hi <- hist(cs, -1:max(cs), plot=FALSE)$intensities
  if (!cumulative) {
    res <- hi
  } else {
    res <- rev(cumsum(rev(hi)))
  }
  
  res
}

closeness <- function(graph, v=igraph.vs.all(graph), mode="all") {
  if (is.character(mode)) {
    mode <- switch(mode, "out"=1, "in"=2, "all"=3)
  }
  
  .Call("R_igraph_closeness", graph, as.igraph.vs(graph, v), as.numeric(mode),
        PACKAGE="igraph")
}

shortest.paths <- function(graph, v=igraph.vs.all(graph), mode="all") {
  if (is.character(mode)) {
    mode <- switch(mode, "out"=1, "in"=2, "all"=3)
  }

  .Call("R_igraph_shortest_paths", graph, as.igraph.vs(graph, v),
        as.numeric(mode),
        PACKAGE="igraph")
}

get.shortest.paths <- function(graph, from, mode="all") {
  if (is.character(mode)) {
    mode <- switch(mode, "out"=1, "in"=2, "all"=3)
  }

  .Call("R_igraph_get_shortest_paths", graph,
        as.numeric(from), as.numeric(mode),
        PACKAGE="igraph")
}

get.all.shortest.paths <- function(graph, from, mode="all") {
  if (is.character(mode)) {
    mode <- switch(mode, "out"=1, "in"=2, "all"=3)
  }

  .Call("R_igraph_get_all_shortest_paths", graph,
        as.numeric(from), as.numeric(mode),
        PACKAGE="igraph")
}

subcomponent <- function(graph, v, mode="all") {
  if (is.character(mode)) {
    mode <- switch(mode, "out"=1, "in"=2, "all"=3)
  }

  .Call("R_igraph_subcomponent", graph, as.numeric(v), as.numeric(mode),
        PACKAGE="igraph")
}

subgraph <- function(graph, v) {
  .Call("R_igraph_subgraph", graph, as.igraph.vs(graph, v),
        PACKAGE="igraph")
}

simplify <- function(graph, remove.multiple=TRUE,
                     remove.loops=TRUE) {
  .Call("R_igraph_simplify", graph, as.logical(remove.multiple),
        as.logical(remove.loops), PACKAGE="igraph")

##   res <- graph
##   vc <- vcount(res)
##   if (remove.loops && vc > 0) {
##     remove <- numeric()
##     for (i in 0:(vc-1)) {
##       neis <- neighbors(graph, i, "out")
##       loops <- sum(neis==i)
##       if (is.directed(graph)) { loops <- loops*2 }
##       remove <- c(remove, rep(i, loops))
##     }
##     res <- delete.edges(res, remove)  
##   }
##   if (remove.multiple) {
##     remove <- numeric()
##     for (i in 0:(vc-1)) {
##       neis <- neighbors(graph, i, "out")
##       dup <- neis[ duplicated(neis) & neis > i ]
##       l <- sum(neis==i)
##       if (l>2) { dup <- c(dup, rep(i, l/4)) }
##       remove <- c(remove, as.numeric(t(matrix(c(rep(i,length(dup)),
##                                                 dup), nc=2))))
##     }
##     res <- delete.edges(res, remove)  
##   }  
  
##   res
}

betweenness <- function(graph, v=igraph.vs.all(graph), directed=TRUE) {
  
  .Call("R_igraph_betweenness", graph, as.igraph.vs(graph, v),
        as.logical(directed),
        PACKAGE="igraph")
}

edge.betweenness <- function(graph, e=igraph.es.all(graph), directed=TRUE) {

  .Call("R_igraph_edge_betweenness", graph, as.logical(directed),
        PACKAGE="igraph")[ as.vector(e)+1 ]  
}

transitivity <- function(graph, type="undirected") {
  if (is.character(type)) {
    type <- switch(type, "undirected"=0)
  }

  .Call("R_igraph_transitivity", graph, as.numeric(type),
        PACKAGE="igraph")
}
