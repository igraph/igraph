
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

graph.empty.adjacencylist.default <- function(..., type="adjacencylist") {
  
  gal <- list(type="adjacencylist", ...)
  if (!"directed" %in% names(gal)) {
    gal$directed <- TRUE
  }
  realn <- ifelse(is.null(gal$n), 0, gal$n)
  gal$n <- 0
  gal$named <- FALSE
  gal$id <- igraph.i.create.id()
  
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
  res$gal$id <- igraph.i.create.id()
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
  res$gal$id <- igraph.i.create.id()
  
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

delete.edges.adjacencylist.default <- function(graph, edges) {

  delete.edges.common(graph, edges)
  
  res <- graph
  res$gal$id <- igraph.i.create.id()
  
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
  res$gal$id <- igraph.i.create.id()
  v <- as.numeric(v)
  res$data$out <- .Call("REST_neighborlist_delete_vertices", res$data$out, v,
                        graph$gal$named, PACKAGE="igraph")
  if (is.directed(graph)) {
    res$data$inc <- .Call("REST_neighborlist_delete_vertices", res$data$inc, v,
                          graph$gal$named, PACKAGE="igraph")
  }
  res$gal$n <- length(res$data$out)
  if (length(res$val) != 0 && length(v)>0 ) {
    res$val <- res$val[-v]
  }
  
  res
}

###################################################################
# Structure query
###################################################################

vcount.adjacencylist.default <- function(graph) {
  res <- as.double(length(graph$data$out))
  res
}

ecount.adjacencylist.default <- function(graph) {
  res <- .Call("REST_i_adjacencylist_ecount", igraph.c.interface,
               graph, PACAKGE="igraph")

  if (!is.directed(graph)) {
    res <- res/2
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

add.graph.attribute.adjacencylist.default <- function(graph,
                                                      attrname, default=NA) {
  if (!is.null(graph$gal[[attrname]])) {
    stop("attribute '", attrname, "' already present")
  }
  
  res <- graph
  res$gal[[attrname]] <- default
  res
}

delete.graph.attribute.adjacencylist.default <- function(graph, attrname) {
  if (is.null(graph$gal[[attrname]])) {
    stop("no such attribute: ", attrname)
  }

  res <- graph
  res$gal[[attrname]] <- NULL
  res
}

get.graph.attribute.adjacencylist.default <- function(graph, attrname=NULL) {
  if(is.null(attrname)) {
    res <- names(graph$gal)
  } else {
    res <- graph$gal[[attrname]]
  }

  res
}

set.graph.attribute.adjacencylist.default <- function(graph, attrname, value) {
  if (is.null(graph$gal[[attrname]])) {
    stop("no such attribute: ", attrname)
  }
  
  res <- graph
  res$gal[[attrname]] <- value
  res
}

add.vertex.attribute.adjacencylist.default <- function(graph,
                                                       attrname, default=NA) {
  if (!is.null(graph$val[[attrname]])) {
    stop("attribute '", attrname, "' already present")
  }

  res <- graph
  res$val[[attrname]] <- replicate(vcount(graph), default, simplify=FALSE)
  
  res
}

delete.vertex.attribute.adjacencylist.default <- function(graph, attrname) {
  if (is.null(graph$val[[attrname]])) {
    stop("no such attribute: ", attrname)
  }

  res <- graph
  res$val[[attrname]] <- NULL

  res
}

get.vertex.attribute.adjacencylist.default <- function(graph, attrname=NULL,
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

set.vertex.attribute.adjacencylist.default <- function(graph, attrname,
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

add.edge.attribute.adjacencylist.default <- function(graph,
                                                     attrname, default=NA) {
  res <- graph
  
  if (is.null(res$eal)) {
    res$eal <- list()
  }

  if (is.null(res$gal$named) || !res$gal$named) {
    res <- name.edges.adjacencylist.default(res)
  }
  
  res$eal[[attrname]] <- new.env(hash=TRUE)  
  attributes(res$eal[[attrname]])$default <- default

  res
}

delete.edge.attribute.adjacencylist.default <- function(graph, attrname) {
  if (is.null(graph$eal[[attrname]])) {
    stop("No such edge attribute")
  }

  res <- graph
  res$eal[[attrname]] <- NULL
  
  res
}

get.edge.attribute.adjacencylist.default <- function(graph, attrname=NULL,
                                                     from=NULL,
                                                     to=NULL) {
  if (is.null(attrname)) {
    res <- names(graph$eal)
  } else {  
    ind <- as.character(get.edge.names.adjacencylist.default(graph,
                                                             from=from, to=to))
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

set.edge.attribute.adjacencylist.default <- function(graph, attrname,
                                                     from=NULL, to=NULL,
                                                     value) {
  if (is.null(graph$eal[[attrname]])) {
    stop("No such attribute", attrname)
  }
  
  res <- graph
  
  if (is.null(from) && is.null(to)) {
    res$eal[[attrname]] <- new.env(hash=TRUE)
    attributes(res$eal[[attrname]])$default <- value
  } else {
    ind <- as.character(get.edge.names.adjacencylist.default(res,
                                                             from=from, to=to))
    for (i in seq(along=ind)) {
      res$eal[[attrname]][[ ind[i] ]] <- value
    }
  }
  res
}

###################################################################
# Iterators
###################################################################

igraph.iterator.adjacencylist.default <- function(graph, type="vid") {
  if (type=="vid") {
    attr <- list(graph.type=igraph.type(graph),
                 type=type, id=g.a(graph, "id"))
    data <- 1
  } else if (type=="eid") {
    attr <- list(graph.type=igraph.type(graph),
                 type=type, id=g.a(graph, "id"))
    ec <- ecount(graph)
    if (ec==0) {
      data <- c(-1, -1)
    } else {
      vid <- 1
      nei <- 1
      if (is.directed(graph)) {
        while (length(graph$data$out[[vid]])==0) {
          vid <- vid + 1
        }
      } else {
        neis <- graph$data$out[[vid]]
        nei <- which(neis >= vid)[1]
        while (is.na(nei)) {
          vid <- vid + 1
          neis <- graph$data$out[[vid]]
          nei <- which(neis >= vid)
        }
      }
      data <- c(vid, nei)
    }
  } else {
    stop("Unknown iterator type")
  }
  
  res <- list(data=data, attr=attr) 
  class(res) <- "igraph.iterator"
  
  res
}

igraph.next.adjacencylist.default <- function(graph, it) {
  if (g.a(graph, "id")!=it$attr$id) {
    stop("Invalid iterator")
  }
  res <- it
  if (res$attr$type == "vid") {
    res$data <- res$data+1
    if (res$data > vcount(graph)+1) {
      stop("No more elements")
    }
  } else if (res$attr$type == "eid") {
    if (res$data[1] == -1) {
      stop("No more elements")
    }
    if (is.directed(graph)) {
      res$data[2] <- res$data[2] + 1
      while (res$data[1] <= vcount(graph) &&
             length(graph$data$out[[res$data[1]]]) < res$data[2]) {
        res$data[1] <- res$data[1] + 1
        res$data[2] <- 1
      }
      if (res$data[1] > vcount(graph)) {
        res$data <- c(-1, -1)
      }
    } else {
      vid <- res$data[1]
      nei <- res$data[2]+1
      l <- FALSE
      while (!l && vid <= vcount(graph)) {
        neis <- graph$data$out[[vid]]
        if (length(neis) < nei) {
          vid <- vid + 1
          nei <- 1
        } else if (length(neis) >= nei && neis[nei] > vid) {
          l <- TRUE
        } else if (neis[nei] < vid) {
          nei <- nei + 1
        } else if (neis[nei] == vid &&
                   which(which(neis == vid) == nei) %% 2 == 0) {
          l <- TRUE
        } else {
          nei <- nei + 1
        }
      }
      if (vid > vcount(graph)) {
        res$data <- c(-1, -1)
      } else {
        res$data <- c(vid, nei)
      }
    }
  } else {
    stop("Invalid iterator type")
  }
  
  res
}

igraph.prev.adjacencylist.default <- function(graph, it) {
  if (g.a(graph, "id")!=it$attr$id) {
    stop("Invalid iterator")
  }
  res <- it
  if (it$attr$type == "vid") {
    res$data <- res$data-1
    if (res$data==0) {
      stop("No previous element")
    } 
#  } else if (it$attr$type == "eid") {
    # TODO
  } else {
    stop("Invalid iterator type")
  }
  
  res
}

igraph.end.adjacencylist.default <- function(graph, it) {
  if (g.a(graph, "id") != it$att$id) {
    stop("Invalid iterator")
  }
  if (it$attr$type == "vid") {
    res <- (it$data==vcount(graph)+1)
  } else if (it$attr$type == "eid") {
    res <- (it$data[1] == -1)
  } else {
    stop("Invalid iterator type")
  }

  res
}

igraph.get.adjacencylist.default <- function(graph, it) {
  if (g.a(graph, "id") != it$att$id) {
    stop("Invalid iterator")
  }
  if (igraph.end(graph, it)) {
    stop("No more elements")
  }
  if (it$attr$type == "vid") {
    res <- it$data
  } else if (it$attr$type == "eid") {
    res <- unname(c(it$data[1], graph$data$out[[it$data[1]]][it$data[2]]))
  } else {
    stop("Invalid iterator type")
  }
  
  res
} 

igraph.getattr.adjacencylist.default <- function(graph, it, attr=NULL) {
  if (g.a(graph, "id") != it$att$id) {
    stop("Invalid iterator")
  }

  if (igraph.end(graph, it)) {
    stop("No more elements")
  }
  if (it$attr$type == "vid") {
    if (is.null(attr)) {
      res <- lapply(graph$val, "[[", it$data)
    } else {
      res <- graph$val[[attr]][[it$data]]
    }
  } else if (it$attr$type == "eid") {
    if (is.null(attr)) {
      res <- list()
      for (n in names(graph$eal)) {
        res[[n]] <-
          graph$eal[[n]][[names(graph$data$out[[it$data[1]]])[it$data[2]]]]
        if (is.null(res[[n]])) {
          res[[n]] <- attributes(graph$eal[[n]])$default
        }
      }
    } else {
      res <-
        graph$eal[[attr]][[names(graph$data$out[[it$data[1]]])[it$data[2]]]]
    }
  } else {
    stop("Invalid iterator type")
  }
  
  res
}

###################################################################
# Reimplementations for speedup
###################################################################

degree.adjacencylist <- function(graph, v=1:vcount(graph), mode="total",
                                 loops=TRUE) {

  res <- numeric(length(v))
  if (!is.directed(graph) || mode %in% c("total", "out")) {
    if (loops) {
      res <- res + sapply(graph$data$out[v], length)
    } else {
      res <- res + sapply(sapply(seq(along=v), function(a)
                                 graph$data$out[[a]][graph$data$out[[a]]!=a]),
                          length)
    }
  }
  if (is.directed(graph) && mode %in% c("total", "in")) {
    if (loops) {
      res <- res + sapply(graph$data$inc[v], length)
    } else {
      res <- res + sapply(sapply(seq(along=v), function(a)
                                 graph$data$inc[[a]][graph$data$inc[[a]]!=a]),
                          length)
    }
  }
  
  res
}

###################################################################
# Internal functions
###################################################################

remove.one <- function(v, e, silent=FALSE) {
  idx <- (1:length(v))[v==e][1]
  if (is.na(idx) && !silent) {
    stop("No such element to remove")
  }
  v[-idx]
}

name.edges.adjacencylist.default <- function(graph) {

  res <- graph

  ind <- 1
  if (is.directed(res)) { # Directed
    for (i in 1:length(res$data$inc)) {
      names(res$data$inc[[i]]) <- numeric()
    }
    for (i in 1:length(res$data$out)) {
      deg <- length(res$data$out[[i]])
      names(res$data$out[[i]]) <- seq(from=ind, length=deg)

      for (j in seq(along=res$data$out[[i]])) {
        to <- res$data$out[[i]] [j]
        ix <- (1:length(res$data$inc[[to]]))[ res$data$inc[[to]] == i &
                                             is.na(names(res$data$inc[[to]])) ]
        names(res$data$inc[[to]])[ix] <- ind
        ind <- ind + 1
      }
    }
    
  } else { # Undirected
    for (i in 1:length(res$data$out)) {
      names(res$data$out[[i]]) <- numeric()
    }
    for (i in 1:length(res$data$out)) {
      for (j in seq(along=res$data$out[[i]])) {
        if (!is.na(names(res$data$out[[i]])[j])) { next }
        to <- res$data$out[[i]] [j]
        names(res$data$out[[i]])[j] <- ind
        ix <- (1:length(res$data$out[[to]]))[ res$data$out[[to]] == i &
                                             is.na(names(res$data$out[[to]])) ]
        names(res$data$out[[to]])[ix] <- ind
        ind <- ind + 1
      }
    }
  }
  
  res$gal$named <- TRUE
  attributes(res$data)$lastname <- ind-1
  
  res
}

get.edge.names.adjacencylist.default <- function(graph, from=NULL, to=NULL) {

  if (is.null(from) && is.null(to)) {
    res <- unlist(lapply(graph$data$out, names))
    if (!is.directed(graph)) res <- unique(res)
  } else if (is.null(to)) {
    res <- names(graph$data$out[[from]])
    if (!is.directed(graph)) res <- unique(res)
  } else if (is.null(from)) {
    if (is.directed(graph)) {
      res <- names(graph$data$inc[[to]])
    } else {
      res <- unique(names(graph$data$out[[to]]))
    }
  } else {
    res <- unlist(lapply(graph$data$out[from], function(v) {
      names(v) [ v %in% to ]
    }))
    if (!is.directed(graph)) res <- unique(res)
  }

  res
}
