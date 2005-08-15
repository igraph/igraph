
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
# Building graphs 
###################################################################

graph.empty.edgelist.default <- function(..., type="edgelist") {

  gal <- list(type="edgelist", ...)
  if (!"directed" %in% names(gal)) {
    gal$directed <- TRUE
  }
  realn <- ifelse(is.null(gal$n), 0, gal$n)
  gal$n <- 0

  data <- list(data=matrix(0, nc=2, nr=0))
  val <- list()
  eal <- list()
  
  res <- list(data=data, gal=gal, val=val, eal=eal)
  class(res) <- "graph"
  res <- add.vertices(res, realn)

  res
}

add.edges.edgelist.default <- function(graph, edges) {

  add.edges.common(graph, edges)

  res <- graph
  if (length(edges)!=0) {
    edges <- matrix(edges, nc=2, byrow=TRUE)
    if (!is.directed(graph)) {
      edges <- edgelist.switch(edges)
    }    
    res$data$data <- matrix( c(t(res$data$data), t(edges)), nc=2, byrow=TRUE)
    for (n in names(res$eal)) {
      tmp <- attributes(res$eal[[n]])
      res$eal[[n]] <- append(res$eal[[n]], vector(mode="list",
                                         length(edges)/2))
      attributes(res$eal[[n]]) <- tmp
    }
    order <- edgelist.order(res$data$data)
    res$data$data <- res$data$data[order,]
    for (n in names(res$eal)) {
      tmp <- attributes(res$eal[[n]])
      res$eal[[n]] <- res$eal[[n]] [order]
      attributes(res$eal[[n]]) <- tmp
    }
  }

  res
}

add.vertices.edgelist.default <- function(graph, nv) {

  add.vertices.common(graph, nv)

  res <- graph
  res$gal$n <- res$gal$n + nv

  # val
  if (length(res$val) != 0) {
    res$val <- append(res$val, rep(list(list()), nv))
  }
  
  res
}

delete.edges.edgelist.default <- function(graph, edges) {

  delete.edges.common(graph, edges)

  res <- graph
  if (length(edges)>0) {
    edges <- matrix(edges, nc=2, byrow=TRUE)
    if (!is.directed(graph)) {
      edges <- edgelist.switch(edges)
    }    
    edges <- edgelist.sort(edges)
    remaining <- .Call("REST_edgelist_delete_edges_sorted",
                       res$data$data, edges, PACKAGE="igraph")
    res$data$data <- res$data$data[remaining,]
    for (n in names(res$eal)) {
      res$eal[[n]] <- res$eal[[n]][remaining]
    }
  }
  
  res
}

delete.vertices.edgelist.default <- function(graph, v) {

  res <- graph
  v <- as.numeric(v)
  res$data$data <- .Call("REST_edgelist_delete_vertices",
                         res$data$data, res$gal$n, v,
                         PACKAGE="igraph")
  res$gal$n <- res$data$data[[1]]

  for (n in names(res$eal)) {
    tmp <- attributes(res$eal[[n]])
    res$eal[[n]] <- res$eal[[n]] [ res$data$data[[3]] ]
    attributes(res$eal[[n]]) <- tmp
  }
  
  res$data$data <- res$data$data[[2]]

  for (n in names(res$val)) {
    res$val[[n]] <- res$val[[n]] [-v]
  }
  
  res
}

###################################################################
# Structure query
###################################################################

vcount.edgelist.default <- function(graph) {
  res <- graph$gal$n
  res
}

ecount.edgelist.default <- function(graph) {
  res <- length(graph$data$data)/2
  res
}

neighbors.edgelist.default <- function(graph, v, mode="out") {

  neighbors.common(graph, v, mode)

  v <- as.numeric(v)
  if (is.directed(graph)) {
    mode=switch(mode, "out"=1, "in"=2, "all"=3)
  } else {
    mode=3
  }

  res <- .Call("REST_edgelist_neighbors", graph$data$data, v, as.numeric(mode),
               PACKAGE="igraph")
  
  res
}

###################################################################
# Attributes
###################################################################

add.graph.attribute.edgelist.default <- function(graph,
                                                      attrname, default=NA) {
  if (!is.null(graph$gal[[attrname]])) {
    stop("attribute '", attrname, "' already present")
  }
  
  res <- graph
  res$gal[[attrname]] <- default
  res
}

delete.graph.attribute.edgelist.default <- function(graph, attrname) {
  if (is.null(graph$gal[[attrname]])) {
    stop("no such attribute: ", attrname)
  }

  res <- graph
  res$gal[[attrname]] <- NULL
  res
}

get.graph.attribute.edgelist.default <- function(graph, attrname=NULL) {
  if(is.null(attrname)) {
    res <- names(graph$gal)
  } else {
    res <- graph$gal[[attrname]]
  }

  res
}

set.graph.attribute.edgelist.default <- function(graph, attrname, value) {
  if (is.null(graph$gal[[attrname]])) {
    stop("no such attribute: ", attrname)
  }
  
  res <- graph
  res$gal[[attrname]] <- value
  res
}

add.vertex.attribute.edgelist.default <- function(graph,
                                                       attrname, default=NA) {
  if (!is.null(graph$val[[attrname]])) {
    stop("attribute '", attrname, "' already present")
  }

  res <- graph
  res$val[[attrname]] <- replicate(vcount(graph), default, simplify=FALSE)
  
  res
}

delete.vertex.attribute.edgelist.default <- function(graph, attrname) {
  if (is.null(graph$val[[attrname]])) {
    stop("no such attribute: ", attrname)
  }

  res <- graph
  res$val[[attrname]] <- NULL

  res
}

get.vertex.attribute.edgelist.default <- function(graph, attrname=NULL,
                                                       v=NULL) {
  if (is.null(attrname)) {
    res <- names(graph$val)
  } else {
    if (is.null(graph$val[[attrname]])) {
      stop("no such attribute: ", attrname)
    }
    
    if (is.null(v)) {
      res <- graph$val[[attrname]]
    } else {
      res <- graph$val[[attrname]][v]
    }
    if (length(res)==1) {
      res <- res[[1]]
    }
  }

  res
}

set.vertex.attribute.edgelist.default <- function(graph, attrname,
                                                  v=NULL, value) {
  if (is.null(graph$val[[attrname]])) {
    stop("no such attribute: ", attrname)
  }

  res <- graph
  if (is.null(v)) {
    v <- 1:vcount(graph)
  }
  res$val[[attrname]][v] <- value
  
  res
}

add.edge.attribute.edgelist.default <- function(graph,
                                                attrname, default=NA) {
  res <- graph

  if (is.null(res$eal)) {
    res$eal <- list()
  }

  res$eal[[attrname]] <- vector(mode="list", ecount(graph))
  attributes(res$eal[[attrname]])$default <- default

  res
}

delete.edge.attribute.edgelist.default <- function(graph, attrname) {
  if (is.null(graph$eal[[attrname]])) {
    stop("No such edge attribute")
  }

  res <- graph
  res$eal[[attrname]] <- NULL

  res
}

get.edge.attribute.edgelist.default <- function(graph, attrname=NULL,
                                                from=NULL,
                                                to=NULL) {
  if (is.null(attrname)) {
    res <- names(graph$eal)
  } else {
    ind <- edgelist.get.edge.indices(graph, from=from, to=to)
    res <- vector(mode="list", length(ind))
    default <- attributes(graph$eal[[attrname]])$default
    for (i in seq(along=ind)) {
      tmp <- graph$eal[[attrname]][[ ind[i] ]]
      res[[i]] <- ifelse(is.null(tmp), default, tmp)
    }

    if (length(res)==1) {
      res <- res[[1]]
    }
  }

  res
}

set.edge.attribute.edgelist.default <- function(graph, attrname,
                                                from=NULL, to=NULL,
                                                value) {
  if (is.null(graph$eal[[attrname]])) {
    stop("No such attribute", attrname)
  }

  res <- graph

  if (is.null(from) && is.null(to)) {
    res$eal[[attrname]] <- vector(mode="list", ecount(graph))
    attributes(res$eal[[attrname]])$default <- value
  } else {
    ind <- edgelist.get.edge.indices(res, from=from, to=to)
    for (i in seq(along=ind)) {
      res$eal[[attrname]][[ ind[i] ]] <- value
    }
  }
  res
}

###################################################################
# Internal
###################################################################

edgelist.switch <- function(edges) {
  chidx <- edges[,1] > edges[,2]
  tmp <- edges[chidx,1]
  edges[chidx,1] <- edges[chidx,2]
  edges[chidx,2] <- tmp
  edges
}

edgelist.order <- function(edges) {
  order <- order(edges[,1], edges[,2])
  order
}

edgelist.sort <- function(edges) {
  edges <- t(rbind(edges[,1], edges[,2]) [, order(edges[,1], edges[,2]) ])
  edges
}

edgelist.get.edge.indices <- function(graph, from=NULL, to=NULL) {

  if (is.null(from) && is.null(to)) {
    res <- seq(along=graph$eal[[1]])
  } else if (is.null(to)) {
    if (is.directed(graph)) {
      res <- which(graph$data$data[,1] %in% from)
    } else {
      res <- which(graph$data$data[,1] %in% from |
                   graph$data$data[,2] %in% from)
    }
  } else if (is.null(from)) {
    if (is.directed(graph)) {
      res <- which(graph$data$data[,2] %in% to)
    } else {
      res <- which(graph$data$data[,1] %in% to |
                   graph$data$data[,2] %in% to)
    }
  } else {
    if (is.directed(graph)) {
      res <- which(graph$data$data[,1] %in% from &
                   graph$data$data[,2] %in% to)
    } else {
      res <- which((graph$data$data[,1] %in% from &
                    graph$data$data[,2] %in% to) |
                   (graph$data$data[,1] %in% to &
                    graph$data$data[,2] %in% from))
    }
  }

  res
}
