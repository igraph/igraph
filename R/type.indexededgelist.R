
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

# TODO: directed & undirected case

###################################################################
# Building graphs 
###################################################################

graph.empty.indexededgelist.default <- function(...,
                                                type="indexededgelist") {
  gal <- list(type="indexededgelist", ...)
  if (!"directed" %in% names(gal)) {
    gal$directed <- TRUE
  }
  realn <- ifelse(is.null(gal$n), 0, gal$n)
  gal$n <- 0
  gal$id <- igraph.i.create.id()

  data <- list(el=matrix(0, nc=2, nr=0),   # edge list matrix
               oi=numeric(0),              # out index (unneeded if sorted)
               ii=numeric(0),              # in index
               os=1,                       # out start index
               is=1)                       # in start index
  
  val <- list()
  eal <- list()

  res <- list(data=data, gal=gal, val=val, eal=eal)
  class(res) <- "graph"
  res <- add.vertices(res, realn)

  res
}

add.edges.indexededgelist.default <- function(graph, edges) {

  add.edges.common(graph, edges)

  res <- graph
  res$gal$id <- igraph.i.create.id()
  if (length(edges)!=0) {
    edges <- matrix(edges, nc=2, byrow=TRUE)
    res$data$el <- matrix( c(t(res$data$el), t(edges)), nc=2,
                          byrow=TRUE)
     for (n in names(res$eal)) {
       tmp <- attributes(res$eal[[n]])
       res$eal[[n]] <- append(res$eal[[n]], vector(mode="list",
                                                   length(edges)/2))
       attributes(res$eal[[n]]) <- tmp
     }
    res$data$oi <- as.numeric(order(edges[,1], edges[,2]))
    res$data$ii <- as.numeric(order(edges[,2], edges[,1]))
    res$data$os <- indexededgelist.create.out.index(res)
    res$data$is <- indexededgelist.create.in.index(res)
  }
  
  res
}

 add.vertices.indexededgelist.default <- function(graph, nv) {

   add.vertices.common(graph, nv)

   res <- graph
   res$gal$id <- igraph.i.create.id()
   res$gal$n <- res$gal$n + nv

   # val
   if (length(res$val) != 0) {
     res$val <- append(res$val, rep(list(list()), nv))
   }

   # update start indices
   ec <- ecount(graph)
   res$data$os <- c(res$data$os, rep(ec+1, nv))
   res$data$is <- c(res$data$is, rep(ec+1, nv))

   res
 }

 delete.edges.indexededgelist.default <- function(graph, edges) {
   ## TODO
   stop("Not implemented yet")
 }

 delete.vertices.indexededgelist.default <- function(graph, v) {
   ## TODO
   stop("Not implemented yet")
 }

 ###################################################################
 # Structure query
 ###################################################################

 vcount.indexededgelist.default <- function(graph) {
   res <- graph$gal$n
   res
 }

 ecount.indexededgelist.default <- function(graph) {
   res <- length(graph$data$el)/2
   res
 }

 neighbors.indexededgelist.default <- function(graph, v, mode="out") {

   neighbors.common(graph, v, mode)

   v <- as.numeric(v)

   res <- numeric()
   if (mode %in% c("in", "all")) {
     if (graph$data$is[v] < graph$data$is[v+1]) {
       idx <- graph$data$is[v] : (graph$data$is[v+1]-1)
       idx <- graph$data$ii[idx]
       res <- c(res, graph$data$el[idx, 1])
     }
   }
   if (mode %in% c("out", "all")) {
     if (graph$data$os[v] < graph$data$os[v+1]) {
       idx <- graph$data$os[v] : (graph$data$os[v+1]-1)
       idx <- graph$data$oi[idx]
       res <- c(res, graph$data$el[idx, 2])
     }
   }

   res
 }

###################################################################
# Attributes
###################################################################

add.graph.attribute.indexededgelist.default <- function(graph,
                                                        attrname,
                                                        default=NA) { 
  if (!is.null(graph$gal[[attrname]])) {
    stop("attribute '", attrname, "' already present")
  }
  
  res <- graph
  res$gal[[attrname]] <- default
  res
}

delete.graph.attribute.indexededgelist.default <- function(graph,
                                                           attrname) { 
  if (is.null(graph$gal[[attrname]])) {
    stop("no such attribute: ", attrname)
  }

  res <- graph
  res$gal[[attrname]] <- NULL
  res
}

get.graph.attribute.indexededgelist.default <- function(graph,
                                                        attrname=NULL) { 
  if(is.null(attrname)) {
    res <- names(graph$gal)
  } else {
    res <- graph$gal[[attrname]]
  }

  res
}

set.graph.attribute.indexededgelist.default <- function(graph,
                                                        attrname, value) { 
  if (is.null(graph$gal[[attrname]])) {
    stop("no such attribute: ", attrname)
  }
  
  res <- graph
  res$gal[[attrname]] <- value
  res
}

add.vertex.attribute.indexededgelist.default <- function(graph,
                                                         attrname,
                                                         type="simple",
                                                         default=NA) {
  if (!is.null(graph$val[[attrname]])) {
    stop("attribute '", attrname, "' already present")
  }

  res <- graph
  vc <- vcount(graph)
  if (type=="complex") {
    res$val[[attrname]] <- replicate(vc, default, simplify=FALSE)
  } else {
    res$val[[attrname]] <- rep(default, vc)
  }
  
  res
}

delete.vertex.attribute.indexededgelist.default <- function(graph,
                                                            attrname) { 
  if (is.null(graph$val[[attrname]])) {
    stop("no such attribute: ", attrname)
  }

  res <- graph
  res$val[[attrname]] <- NULL

  res
}

get.vertex.attribute.indexededgelist.default <- function(graph, attrname=NULL,
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

set.vertex.attribute.indexededgelist.default <- function(graph, attrname,
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

add.edge.attribute.indexededgelist.default <- function(graph,
                                                       attrname,
                                                       type="complex",
                                                       default=NA) {
  res <- graph

  if (is.null(res$eal)) {
    res$eal <- list()
  }

  ec <- ecount(res)
  if (type=="complex") {
    res$eal[[attrname]] <- replicate(ec, default, simplify=FALSE)
  } else if (type=="simple") {
    res$eal[[attrname]] <- rep(default, ec)
  }
  attributes(res$eal[[attrname]])$default <- default

  res    
}

delete.edge.attribute.indexededgelist.default <- function(graph,
                                                          attrname) {
  if (is.null(graph$eal[[attrname]])) {
    stop("No such edge attribute")
  }

  res <- graph
  res$eal[[attrname]] <- NULL

  res
}

get.edge.attribute.indexededgelist.default <- function(graph,
                                                       attrname=NULL,
                                                       from=NULL,
                                                       to=NULL) {
  if (is.null(attrname)) {
    res <- names(graph$eal)
  } else {
    idx <- indexededgelist.get.edge.index(graph, from, to)
    res <- graph$eal[[attrname]] [ idx ]
  }
  res
}

set.edge.attribute.indexededgelist.default <- function(graph,
                                                       attrname,
                                                       from=NULL,
                                                       to=NULL,
                                                       value) {
  if (is.null(graph$eal[[attrname]])) {
    stop("No such attribute", attrname)
  }
  
  idx <- indexededgelist.get.edge.index(graph, from, to)
  res <- graph
  res$eal[[attrname]] [ idx ] <- value
  
  res
}

###################################################################
# Iterators
###################################################################

igraph.iterator.indexededgelist.default <- function(graph, type="vid") {

  res <- list(type=type, graph.type=igraph.type(graph),
              graph.id=g.a(graph, "id"))
  if (type=="vid") {
    res[[1]] <- 1
    res[[2]] <- igraph.indexededgelist.vid.next
    res[[3]] <- igraph.indexededgelist.vid.end
    res[[4]] <- igraph.indexededgelist.vid.get
    res[[5]] <- igraph.indexededgelist.vid.prev
    res[[6]] <- igraph.indexededgelist.vid.getattr
  } else if (type=="eid") {
    res[[1]] <- 1
    res[[2]] <- igraph.indexededgelist.eid.next
    res[[3]] <- igraph.indexededgelist.eid.end
    res[[4]] <- igraph.indexededgelist.eid.get
    res[[5]] <- igraph.indexededgelist.eid.prev
    res[[6]] <- igraph.indexededgelist.eid.getattr
  }
  
  class(res) <- "igraph.iterator"  
  res
}

igraph.indexededgelist.vid.next <- function(graph, it) {
  it[[1]] <- it[[1]] + 1
  it
}

igraph.indexededgelist.vid.prev <- function(graph, it) {
  it[[1]] <- it[[1]] - 1
  it
}

igraph.indexededgelist.vid.end <- function(graph, it) {
  it[[1]] > graph$gal$n
}

igraph.indexededgelist.vid.get <- function(graph, it) {
  it[[1]]
}

igraph.indexededgelist.vid.getattr <- function(graph, it, attr=NULL) {
  if (is.null(attr)) {
    lapply(graph$val, "[[", it[[1]])
  } else {
    graph$val[[attr]][[it[[1]]]]
  }   
}

igraph.indexededgelist.eid.next <- function(graph, it) {
  it[[1]] <- it[[1]] + 1
  it
}

igraph.indexededgelist.eid.prev <- function(graph, it) {
  it[[1]] <- it[[1]] - 1
  it
}

igraph.indexededgelist.eid.end <- function(graph, it) {
  it[[1]] > nrow(graph$data$el)
}

igraph.indexededgelist.eid.get <- function(graph, it) {
  graph$data$el[ it[[1]], ]
}

igraph.indexededgelist.eid.getattr <- function(graph, it, attr=NULL) {
  if (is.null(attr)) {
    lapply(graph$eal, "[", it[[1]])
  } else {
    graph$eal[[attr]][it[[1]]]
  }   
}

###################################################################
# Internal functions
###################################################################

indexededgelist.create.out.index <- function(graph) {

  res <- .Call("REST_indexededgelist_create_startindex",
               graph$data$el[,1], graph$data$oi, graph$gal$n,
               PACKAGE="igraph")
  res
}

indexededgelist.create.in.index <- function(graph) {

  res <- .Call("REST_indexededgelist_create_startindex",
               graph$data$el[,2], graph$data$ii, graph$gal$n,
               PACKAGE="igraph")
  res
}

indexededgelist.get.edge.index <- function(graph, from=NULL, to=NULL)
{
  ## TODO: more edges at the same time
  if (is.null(from) && is.null(to)) {
    res <- seq(along=graph$data$el[,1])
  } else if(is.null(from)) {
    if(graph$data$is[to] < graph$data$is[to+1]) {
      res <- graph$data$ii[ graph$data$is[to] :
                           (graph$data$is[to+1]-1) ]
    } else {
      res <- numeric()
    }
  } else if (is.null(to)) {
    if (graph$data$os[from] < graph$data$os[from+1]) {
      res <- graph$data$oi[ graph$data$os[from] :
                           (graph$data$os[from+1]-1) ]
    } else {
      res <- numeric()
      }
  } else {
    idx <- graph$data$os[from] : (graph$data$os[from+1]-1)
    idx <- idx [graph$data$el[ graph$data$oi[idx], 2] == to]
    res <- graph$data$oi[idx]
  }
  
  res
}
