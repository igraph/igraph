
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

graph.empty.adjacencylist.default <- function(type="adjacencylist", ...) {
  
  gal <- list(type="adjacencylist", ...)
  if (!"directed" %in% names(gal)) {
    gal$directed <- TRUE
  }
  realn <- ifelse(is.null(gal$n), 0, gal$n)
  gal$n <- 0
  gal$named <- FALSE
  
  data <- list(out=list(),inc=list())
  val <- list()
  eal <- list()
  
  res <- list(data=data, gal=gal, val=val, eal=eal)
  class(res) <- "graph"
  res <- add.vertices(res, realn)
  
  res
}

add.edges.adjacencylist.default <- function(graph, edges) {

  add.edges.common(graph, edges)

  res <- graph
  if (length(edges) != 0) {
    res$data <- .Call("REST_add_edges_adjacencylist", res$data,
                      as.numeric(edges), is.directed(graph), res$gal$named,
                      attributes(res$data)$lastname,
                      PACKAGE="igraph")
  }
  
  res
}

add.vertices.adjacencylist.default <- function(graph, nv) {

  add.vertices.common(graph, nv)

  res <- graph

  # data
  res$gal$n <- res$gal$n + nv
  res$data$out <- c(res$data$out, replicate(nv, numeric()))
  if (is.directed(graph)) {
    res$data$inc <- c(res$data$inc, replicate(nv, numeric()))
  }

  # val
  if (length(res$val) != 0) {
    res$val <- append(res$val, rep(list(list()), nv))
  }

  # keep names up to date
  if (res$gal$named && nv > 0) {
    for (i in 1:nv) {
      names(res$data$out[[graph$gal$n+i]]) <- numeric()
      if (is.directed(graph)) {
        names(res$data$inc[[graph$gal$n+i]]) <- numeric()
      }
    }
  }
  
  res
}

remove.one <- function(v, e, silent=FALSE) {
  idx <- (1:length(v))[v==e][1]
  if (is.na(idx) && !silent) {
    stop("No such element to remove")
  }
  v[-idx]
}

delete.edges.adjacencylist.default <- function(graph, edges) {

  delete.edges.common(graph, edges)
  
  res <- graph
  
  if (length(edges) > 0) {
    for (i in seq(1, length(edges), by=2)) {
      res$data$out[[ edges[i] ]] <- remove.one(res$data$out[[ edges[i] ]],
                                               edges[i+1])
      if (is.directed(res)) {
        res$data$inc[[ edges[i+1] ]] <- remove.one(res$data$inc[[edges[i+1] ]],
                                                       edges[i])
      } else {
        res$data$out[[ edges[i+1] ]] <- remove.one(res$data$out[[edges[i+1] ]],
                                                       edges[i])
      }
    }
  }
  
  res  
}

delete.vertices.adjacencylist.default <- function(graph, v) {

  res <- graph
  v <- as.numeric(v)
  res$data$out <- .Call("REST_neighborlist_delete_vertices", res$data$out, v,
                        graph$gal$named, PACKAGE="simplegraph")
  if (is.directed(graph)) {
    res$data$inc <- .Call("REST_neighborlist_delete_vertices", res$data$inc, v,
                          graph$gal$named, PACKAGE="simplegraph")
  }
  res$gal$n <- length(res$data$out)
  if (length(res$val) != 0 && length(v)>0 ) {
    res$val <- res$val[-v]
  }
  
  if (graph.verbose) {
    print("debug: `delete.vertices.adjacencylist.default' finished")
  }
  res
}

###################################################################
# Structure query
###################################################################

vcount.adjacencylist.default <- function(graph) {
  res <- length(graph$data$out)
  res
}

ecount.adjacencylist.default <- function(graph) {
  if (length(graph)==0) {
    res <- 0
  } else {
    res <- sum(sapply(graph$data$out, length))
    if (!is.directed(graph)) {
      res <- res/2
    }
  }
  res
}

neighbors.adjacencylist.default <- function(graph, v, mode="out") {

  neighbors.common(graph, v, mode)

  if (is.directed(graph)) {
    res <- numeric()
    if (mode %in% c("out", "all")) {
      res <- c(res, graph$data$out[[v]])
    }
    if (mode %in% c("in", "all")) {
      res <- c(res, graph$data$inc[[v]])
    } 
  } else {
    # undirected
    res <- graph$data$out[[v]]
  }
  res
}

###################################################################
# Attributes
###################################################################

# TODO
add.graph.attribute.adjacencylist.default <- function(graph,
                                                      attrname, default=NA) {
}

# TODO
delete.graph.attribute.adjacencylist.default <- function(graph, attrname) {
}

get.graph.attribute.adjacencylist.default <- function(graph, attrname) {
  res <- graph$gal[[attrname]]
  res
}

# TODO
set.graph.attribute.adjacencylist.default <- function(graph, attrname, value) {
}

# TODO
add.vertex.attribute.adjacencylist.default <- function(graph,
                                                       attrname, default=NA) {
}

# TODO
delete.vertex.attribute.adjacencylist.default <- function(graph, attrname) {
}

# TODO
get.vertex.attribute.adjacencylist.default <- function(graph, attrname, v) {
}

# TODO
set.vertex.attribute.adjacencylist.default <- function(graph, attrname,
                                                       v, value) {
}

# TODO
add.edge.attribute.adjacencylist.default <- function(graph,
                                                     attrname, default=NA) {
}

# TODO
delete.edge.attribute.adjacencylist.default <- function(graph, attrname) {
}

# TODO
get.edge.attribute.adjacencylist.default <- function(graph, attrname,
                                                     from=NULL, to=NULL) {
}

# TODO
set.edge.attribute.adjacencylist.default <- function(graph, attrname,
                                                     from=NULL, to=NULL,
                                                     value) {
}
