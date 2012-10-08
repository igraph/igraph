
#   IGraph R package
#   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
#   334 Harvard street, Cambridge, MA 02139 USA
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
  if (!all(sapply(graphs, is.igraph))) {
    stop("Not a graph object")
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_disjoint_union", graphs,
               PACKAGE="igraph")

  ## Graph attributes
  galist <- lapply(graphs, list.graph.attributes)
  gan <- unique(unlist(galist))

  attr <- list()
  for (name in gan) { 
    w <- which(sapply(galist, function(x) name %in% x))
    if (length(w)==1) {
      attr[[name]] <- get.graph.attribute(graphs[[w]], name)
    } else {
      for (w2 in w) {
        nname <- paste(name, sep="_", w2)
        attr[[nname]] <- get.graph.attribute(graphs[[w2]], name)
      }      
    }
  }
  graph.attributes(res) <- attr
    
  ## Vertex attributes
  attr <- list()
  vc <- sapply(graphs, vcount)
  cumvc <- c(0, cumsum(vc))
  for (i in seq_along(graphs)) {
    va <- vertex.attributes(graphs[[i]])
    exattr <- intersect(names(va), names(attr))   # existing and present
    noattr <- setdiff(names(attr), names(va))     # existint and missing
    newattr <- setdiff(names(va), names(attr))    # new
    for (a in seq_along(exattr)) {
      attr[[ exattr[a] ]] <- c(attr[[ exattr[a] ]], va[[ exattr[a] ]])
    }
    for (a in seq_along(noattr)) {
      attr[[ noattr[a] ]] <- c(attr[[ noattr[a] ]], rep(NA, vc[i]))
    }
    for (a in seq_along(newattr)) {
      attr[[ newattr[a] ]] <- c(rep(NA, cumvc[i]), va[[ newattr[a] ]])
    }
  }
  vertex.attributes(res) <- attr

  ## Edge attributes
  attr <- list()
  ec <- sapply(graphs, ecount)
  cumec <- c(0, cumsum(ec))
  for (i in seq_along(graphs)) {
    ea <- edge.attributes(graphs[[i]])
    exattr <- intersect(names(ea), names(attr))   # existing and present
    noattr <- setdiff(names(attr), names(ea))     # existint and missing
    newattr <- setdiff(names(ea), names(attr))    # new
    for (a in seq_along(exattr)) {
      attr[[ exattr[a] ]] <- c(attr[[ exattr[a] ]], ea[[ exattr[a] ]])
    }
    for (a in seq_along(noattr)) {
      attr[[ noattr[a] ]] <- c(attr[[ noattr[a] ]], rep(NA, ec[i]))
    }
    for (a in seq_along(newattr)) {
      attr[[ newattr[a] ]] <- c(rep(NA, cumec[i]), ea[[ newattr[a] ]])
    }
  }
  edge.attributes(res) <- attr
  
  res
}

"%du%" <- function(x,y) {
  graph.disjoint.union(x,y)
}

graph.union <- function(...) {

  graphs <- unlist(recursive=FALSE, lapply(list(...), function(l) {
    if (is.igraph(l)) list(l) else l
  } ))
  if (!all(sapply(graphs, is.igraph))) {
    stop("Not a graph object")
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
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
  if (!all(sapply(graphs, is.igraph))) {
    stop("Not a graph object")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
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
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_difference", big, small,
        PACKAGE="igraph")
}
    
"%m%" <- function(x,y) {
  graph.difference(x,y)
}

graph.complementer <- function(graph, loops=FALSE) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_complementer", graph, as.logical(loops),
        PACKAGE="igraph")
}

graph.compose <- function(g1, g2) {

  if (!is.igraph(g1) || !is.igraph(g2)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_compose", g1, g2,
        PACKAGE="igraph")
}

"%c%" <- function(x,y) {
  graph.compose(x,y)
}

edge <- function(...) {
  structure(list(...), class="igraph.edge")
}

edges <- edge

vertex <- function(...) {
  structure(list(...), class="igraph.vertex")
}

vertices <- vertex

path <- function(...) {
  structure(list(...), class="igraph.path")
}

`+.igraph` <- function(e1, e2) {
  if (!is.igraph(e1) && is.igraph(e2)) {
    tmp <- e1
    e1 <- e2
    e2 <- tmp
  }
  if (is.igraph(e2)) {
    ## Disjoint union of graphs
    res <- graph.disjoint.union(e1,e2)

  } else if ("igraph.edge" %in% class(e2)) {
    ## Adding edges, possibly with attributes
    ## Non-named arguments define the edges
    if (is.null(names(e2))) {
      toadd <- unlist(e2, recursive=FALSE)
      attr <- list()
    } else {
      toadd <- unlist(e2[names(e2)==""])
      attr <- e2[names(e2)!=""]
    }
    res <- add.edges(e1, as.igraph.vs(e1, toadd), attr=attr)

  } else if ("igraph.vertex" %in% class(e2)) {
    ## Adding vertices, possibly with attributes
    ## If there is a single unnamed argument, that contains the vertex names
    wn <- which(names(e2)=="")
    if (length(wn)==1) {
      names(e2)[wn] <- "name"
    } else if (is.null(names(e2))) {
    ## No names at all, everything is a vertex name
      e2 <- list(name=unlist(e2, recursive=FALSE))
    } else if (length(wn)==0) {
    ## If there are no non-named arguments, we are fine
    } else {
    ## Otherwise, all unnamed arguments are collected and used as
    ## vertex names
      nn <- unlist(e2[wn], recursive=FALSE)
      e2 <- c(list(name=nn), e2[names(e2)!=""])
    }
    la <- unique(sapply(e2, length))
    res <- add.vertices(e1, la, attr=e2)

  } else if ("igraph.path" %in% class(e2)) {
    ## Adding edges along a path, possibly with attributes
    ## Non-named arguments define the edges
    if (is.null(names(e2))) {
      toadd <- unlist(e2, recursive=FALSE)
      attr <- list()
    } else {
      toadd <- unlist(e2[names(e2)==""])
      attr <- e2[names(e2)!=""]
    }
    toadd <- as.igraph.vs(e1, toadd)
    lt <- length(toadd)
    if (lt >= 2) {
      toadd <- c(toadd[1], rep(toadd[2:(lt-1)], each=2), toadd[lt])
      res <- add.edges(e1, toadd, attr=attr)
    } else {
      res <- e1
    }
    
  } else if (is.numeric(e2) && length(e2)==1) {
    ## Adding some isolate vertices
    res <- add.vertices(e1, e2)

  } else if (is.character(e2)) {
    ## Adding named vertices
    res <- add.vertices(e1, length(e2), name=e2)
    
  } else {
    stop("Cannot add unknown type to igraph graph")
  }
  res
}

`-.igraph` <- function(e1, e2) {
  if (missing(e2)) {
    stop("Non-numeric argument to negation operator")
  }
  if (is.igraph(e2)) {
    res <- graph.difference(e1, e2)
  } else if ("igraph.vertex" %in% class(e2)) {
    res <- delete.vertices(e1, unlist(e2, recursive=FALSE))
  } else if ("igraph.edge" %in% class(e2)) {
    res <- delete.edges(e1, unlist(e2, recursive=FALSE))
  } else if ("igraph.path" %in% class(e2)) {
    todel <- unlist(e2, recursive=FALSE)
    lt <- length(todel)
    if (lt >= 2) {
      todel <- paste(todel[-lt], todel[-1], sep="|")
      res <- delete.edges(e1, todel)
    } else {
      res <- e1
    }
  } else if ("igraph.vs" %in% class(e2)) {
    res <- delete.vertices(e1, e2)
  } else if ("igraph.es" %in% class(e2)) {
    res <- delete.edges(e1, e2)
  } else if (is.numeric(e2) || is.character(e2)) {
    res <- delete.vertices(e1, e2)
  } else {
    stop("Cannot substract unknown type from igraph graph")
  }
  res
}
