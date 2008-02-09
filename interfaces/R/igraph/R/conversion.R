
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

get.adjacency <- function(graph, type=c("both", "upper", "lower"),
                          attr=NULL, names=TRUE,
                          binary=FALSE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  type <- igraph.match.arg(type)
  type <- switch(type, "upper"=0, "lower"=1, "both"=2)
  
  if (is.null(attr)) {    
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
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

get.edgelist <- function(graph, names=TRUE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- matrix(.Call("R_igraph_get_edgelist", graph, TRUE,
                      PACKAGE="igraph"), nc=2)
  if (names && "name" %in% list.vertex.attributes(graph)) {
    res <- matrix(V(graph)$name[ res+1 ], nc=2)
  }

  res
}

as.directed <- function(graph, mode=c("mutual", "arbitrary")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "arbitrary"=0, "mutual"=1)
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_to_directed", graph, as.numeric(mode),
        PACKAGE="igraph")
}

as.undirected <- function(graph, mode=c("collapse", "each")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "each"=0, "collapse"=1)
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_to_undirected", graph, as.numeric(mode),
        PACKAGE="igraph")  
}

get.adjlist <- function(graph, mode=c("all", "out", "in", "total")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  mode <- igraph.match.arg(mode)
  mode <- as.numeric(switch(mode, "out"=1, "in"=2, "all"=3, "total"=3))
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_get_adjlist", graph, mode,
        PACKAGE="igraph")
}

get.adjedgelist <- function(graph, mode=c("all", "out", "in", "total")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  mode <- igraph.match.arg(mode)
  mode <- as.numeric(switch(mode, "out"=1, "in"=2, "all"=3, "total"=3))
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_get_adjedgelist", graph, mode,
        PACKAGE="igraph")
}

