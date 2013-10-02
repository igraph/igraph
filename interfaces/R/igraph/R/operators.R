
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

rename.attr.if.needed <- function(type, graphs, newsize=NULL, maps=NULL,
                                  ignore=character()) {
  listfun <- switch(type, "g"=list.graph.attributes,
                    "v"=list.vertex.attributes, "e"=list.edge.attributes,
                    stop("Internal igraph error"))
  getfun <- switch(type, "g"=get.graph.attribute, "v"=get.vertex.attribute,
                   "e"=get.edge.attribute, stop("Internal igraph error"))
  alist <- lapply(graphs, listfun)
  an <- unique(unlist(alist))
  an <- setdiff(an, ignore)

  getval <- function(which, name) {
    newval <- getfun(graphs[[which]], name)
    if (!is.null(maps)) {
      idx <- rep(NA, newsize)
      midx <- maps[[which]][ maps[[which]] >= 0 ]
      idx[ midx + 1 ] <- seq_along(midx)
      newval <- newval[idx]
    }
    if (!is.null(newsize)) { length(newval) <- newsize }
    newval
  }
  
  attr <- list()
  for (name in an) {
    w <- which(sapply(alist, function(x) name %in% x))
    if (length(w)==1) {
      attr[[name]] <- getval(w, name)
    } else {
      for (w2 in w) {
        nname <- paste(name, sep="_", w2)
        newval <- getval(w2, name)
        attr[[nname]] <-newval
      }      
    }
  }
  attr
}

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
  graph.attributes(res) <- rename.attr.if.needed("g", graphs)
    
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

.igraph.graph.union.or.intersection <- function(call, ..., byname,
                                                keep.all.vertices) {

  graphs <- unlist(recursive=FALSE, lapply(list(...), function(l) {
    if (is.igraph(l)) list(l) else l
  } ))
  if (!all(sapply(graphs, is.igraph))) {
    stop("Not a graph object")
  }
  if (byname != "auto" && !is.logical(byname)) {
    stop("`bynam' must be \"auto\", or logical")
  }
  nonamed <- sum(sapply(graphs, is.named))
  if (byname == "auto") {
    byname <- all(sapply(graphs, is.named))
    if (nonamed != 0 && nonamed != length(graphs)) {
      warning("Some, but not all graphs are named, not using vertex names")
    }
  } else if (byname && nonamed != length(graphs)) {
    stop("Some graphs are not named")
  }

  edgemaps <- length(unlist(lapply(graphs, list.edge.attributes))) != 0
  
  if (byname) {
    allnames <- lapply(graphs, get.vertex.attribute, "name")
    if (keep.all.vertices) {
      uninames <- unique(unlist(allnames))
      newgraphs <- lapply(graphs, function(g) {
        g <- g + setdiff(uninames, V(g)$name)
        permute.vertices(g, match(V(g)$name, uninames))
      })
    } else {
      uninames <- Reduce(intersect, allnames)
      newgraphs <- lapply(graphs, function(g) {
        g <- g - setdiff(V(g)$name, uninames)
        permute.vertices(g, match(V(g)$name, uninames))
      })
    }
    
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    res <- .Call(call, newgraphs, edgemaps, PACKAGE="igraph")
    maps <- res$edgemaps
    res <- res$graph

    ## We might need to rename all attributes
    graph.attributes(res) <- rename.attr.if.needed("g", graphs)
    vertex.attributes(res) <- rename.attr.if.needed("v", graphs,
                                                    vcount(res),
                                                    ignore="name")
    V(res)$name <- uninames

    ## Edges are a bit more difficult, we need a mapping
    if (edgemaps) {
      edge.attributes(res) <- rename.attr.if.needed("e", graphs,
                                                    ecount(res),
                                                    maps=maps)
    }
  } else {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    res <- .Call(call, graphs, edgemaps, PACKAGE="igraph")
    maps <- res$edgemaps
    res <- res$graph

    ## We might need to rename all attributes
    graph.attributes(res) <- rename.attr.if.needed("g", graphs)
    vertex.attributes(res) <- rename.attr.if.needed("v", graphs,
                                                    vcount(res))

    ## Edges are a bit more difficult, we need a mapping
    if (edgemaps) {
      edge.attributes(res) <- rename.attr.if.needed("e", graphs,
                                                    ecount(res),
                                                    maps=maps)
    }
  }

  res
}

graph.union <- function(..., byname="auto") {
  .igraph.graph.union.or.intersection("R_igraph_union", ..., byname=byname,
                                      keep.all.vertices=TRUE)
}

"%u%" <- function(x,y) {
  graph.union(x,y)
}

graph.intersection <- function(..., byname="auto",
                               keep.all.vertices=TRUE) {
  .igraph.graph.union.or.intersection("R_igraph_intersection", ...,
                                      byname=byname,
                                      keep.all.vertices=keep.all.vertices)
}

"%s%" <- function(x,y) {
  graph.intersection(x,y)
}

graph.difference <- function(big, small, byname="auto") {

  if (!is.igraph(big) || !is.igraph(small)) {
    stop("argument is not a graph")
  }
  if (byname != "auto" && !is.logical(byname)) {
    stop("`bynam' must be \"auto\", or logical")
  }
  nonamed <- is.named(big) + is.named(small)
  if (byname == "auto") {
    byname <- nonamed == 2
    if (nonamed == 1) {
      warning("One, but not both graphs are named, not using vertex names")
    }
  } else if (byname && nonamed != 2) {
    stop("Some graphs are not named")
  }
  
  if (byname) {
    bnames <- V(big)$name
    snames <- V(small)$name
    if (any(! snames %in% bnames)) {
      small <- small - setdiff(snames, bnames)
      snames <- V(small)$name
    }
    perm <- match(bnames, snames)
    if (any(is.na(perm))) {
      perm[is.na(perm)] <- seq(from=vcount(small)+1, to=vcount(big))
    }
    big <- permute.vertices(big, perm)

    on.exit(.Call("R_igraph_finalizer", PACKAGE="igraph"))
    res <- .Call("R_igraph_difference", big, small,
                 PACKAGE="igraph")
    permute.vertices(res, match(V(res)$name, bnames))

  } else {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    .Call("R_igraph_difference", big, small,
          PACKAGE="igraph")
  }
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
