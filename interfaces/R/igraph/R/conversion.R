
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
#   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
#   02110-1301 USA
#
###################################################################

get.adjacency <- function(graph, type="both", attr=NULL, names=TRUE,
                          binary=FALSE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (is.character(type)) {
    type <- switch(type, "upper"=0, "lower"=1, "both"=2)
  }
  
  if (is.null(attr)) {    
    res <- .Call("R_igraph_get_adjacency", graph, as.numeric(type),
                 PACKAGE="igraph")
    if (binary) {
      res <- ifelse(res >= 1, 1, 0)
    }
  } else {
    attr <- as.character(attr)
    if (! attr %in% list.edge.attributes(graph)) {
      stop("no such edge attribute")
    }
    res <- matrix(0, nr=vcount(graph), nc=vcount(graph))
    if (is.directed(graph)) {
      for (i in seq(length=ecount(graph))-1) {
        e <- get.edge(graph, i)
        res[ e[1]+1, e[2]+1 ] <- get.edge.attribute(graph, attr, i)
      }
    } else {
      if (type==0) {
        ## upper
        for (i in seq(length=ecount(graph))-1) {
          e <- get.edge(graph, i)
          res[ min(e)+1, max(e)+1 ] <- get.edge.attribute(graph, attr, i)
        }        
      } else if (type==1) {
        ## lower
        for (i in seq(length=ecount(graph))-1) {
          e <- get.edge(graph, i)
          res[ max(e)+1, min(e)+1 ] <- get.edge.attribute(graph, attr, i)
        }        
      } else if (type==2) {
        ## both
        for (i in seq(length=ecount(graph))-1) {
          e <- get.edge(graph, i)
          res[ e[1]+1, e[2]+1 ] <- get.edge.attribute(graph, attr, i)
          if (e[1] != e[2]) {
            res[ e[2]+1, e[1]+1 ] <- get.edge.attribute(graph, attr, i)
          }
        }
      }
    }
  }

  if (names && "name" %in% list.vertex.attributes(graph)) {
    colnames(res) <- rownames(res) <- V(graph)$name
  }
  
  res
}

get.edgelist <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  matrix(.Call("R_igraph_get_edgelist", graph, TRUE,
               PACKAGE="igraph"), nc=2)
}

as.directed <- function(graph, mode="mutual") {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.character(mode)) {
    mode <- switch(mode, "arbitrary"=0, "mutual"=1)
  }
  
  .Call("R_igraph_to_directed", graph, as.numeric(mode),
        PACKAGE="igraph")
}

as.undirected <- function(graph, mode="collapse") {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.character(mode)) {
    mode <- switch(mode, "each"=0, "collapse"=1)
  }
  
  .Call("R_igraph_to_undirected", graph, as.numeric(mode),
        PACKAGE="igraph")  
}
