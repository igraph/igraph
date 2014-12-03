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
                                  maps2=NULL, ignore=character()) {
  listfun <- switch(type, "g"=graph_attr_names,
                    "v"=vertex_attr_names, "e"=edge_attr_names,
                    stop("Internal igraph error"))
  getfun <- switch(type, "g"=graph_attr, "v"=vertex_attr,
                   "e"=edge_attr, stop("Internal igraph error"))
  alist <- lapply(graphs, listfun)
  an <- unique(unlist(alist))
  an <- setdiff(an, ignore)

  getval <- function(which, name) {
    newval <- getfun(graphs[[which]], name)
    if (!is.null(maps)) {
      tmpval <- newval[ maps[[which]] >= 0 ]
      mm <- maps[[which]][ maps[[which]] >= 0 ] + 1
      newval <- rep(NA, newsize)
      newval[mm] <- tmpval
    }
    if (!is.null(maps2)) {
      newval <- newval[ maps2[[which]] + 1 ]
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



#' Disjoint union of graphs
#' 
#' The union of two or more graphs are created. The graphs are assumed to have
#' disjoint vertex sets.
#' 
#' \code{disjoint_union} creates a union of two or more disjoint graphs.
#' Thus first the vertices in the second, third, etc. graphs are relabeled to
#' have completely disjoint graphs. Then a simple union is created. This
#' function can also be used via the \%du\% operator.
#' 
#' \code{graph.disjont.union} handles graph, vertex and edge attributes.  In
#' particular, it merges vertex and edge attributes using the basic \code{c()}
#' function. For graphs that lack some vertex/edge attribute, the corresponding
#' values in the new graph are set to \code{NA}. Graph attributes are simply
#' copied to the result. If this would result a name clash, then they are
#' renamed by adding suffixes: _1, _2, etc.
#' 
#' Note that if both graphs have vertex names (ie. a \code{name} vertex
#' attribute), then the concatenated vertex names might be non-unique in the
#' result. A warning is given if this happens.
#' 
#' An error is generated if some input graphs are directed and others are
#' undirected.
#' 
#' @aliases graph.disjoint.union %du%
#' @param \dots Graph objects or lists of graph objects.
#' @return A new graph object.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @export
#' @keywords graphs
#' @examples
#' 
#' ## A star and a ring
#' g1 <- make_star(10, mode="undirected")
#' V(g1)$name <- letters[1:10]
#' g2 <- make_ring(10)
#' V(g2)$name <- letters[11:20]
#' str(g1 %du% g2)
#' 
disjoint_union <- function(...) {
  
  graphs <- unlist(recursive=FALSE, lapply(list(...), function(l) {
    if (is_igraph(l)) list(l) else l
  } ))
  if (!all(sapply(graphs, is_igraph))) {
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

  if ("name" %in% names(attr) && any(duplicated(attr$name))) {
    warning("Duplicate vertex names in disjoint union")
  }
  
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

#' @export

"%du%" <- function(x,y) {
  disjoint_union(x,y)
}

.igraph.graph.union.or.intersection <- function(call, ..., byname,
                                                keep.all.vertices) {

  graphs <- unlist(recursive=FALSE, lapply(list(...), function(l) {
    if (is_igraph(l)) list(l) else l
  } ))
  if (!all(sapply(graphs, is_igraph))) {
    stop("Not a graph object")
  }
  if (byname != "auto" && !is.logical(byname)) {
    stop("`bynam' must be \"auto\", or logical")
  }
  nonamed <- sum(sapply(graphs, is_named))
  if (byname == "auto") {
    byname <- all(sapply(graphs, is_named))
    if (nonamed != 0 && nonamed != length(graphs)) {
      warning("Some, but not all graphs are named, not using vertex names")
    }
  } else if (byname && nonamed != length(graphs)) {
    stop("Some graphs are not named")
  }

  edgemaps <- length(unlist(lapply(graphs, edge_attr_names))) != 0
  
  if (byname) {
    allnames <- lapply(graphs, vertex_attr, "name")
    if (keep.all.vertices) {
      uninames <- unique(unlist(allnames))
      newgraphs <- lapply(graphs, function(g) {
        g <- g + setdiff(uninames, V(g)$name)
        permute(g, match(V(g)$name, uninames))
      })
    } else {
      uninames <- Reduce(intersect, allnames)
      newgraphs <- lapply(graphs, function(g) {
        g <- g - setdiff(V(g)$name, uninames)
        permute(g, match(V(g)$name, uninames))
      })
    }
    
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    res <- .Call(call, newgraphs, edgemaps, PACKAGE="igraph")
    maps <- res$edgemaps
    res <- res$graph

    ## We might need to rename all attributes
    graph.attributes(res) <- rename.attr.if.needed("g", newgraphs)
    vertex.attributes(res) <- rename.attr.if.needed("v", newgraphs,
                                                    vcount(res),
                                                    ignore="name")
    V(res)$name <- uninames

    ## Edges are a bit more difficult, we need a mapping
    if (edgemaps) {
      edge.attributes(res) <- rename.attr.if.needed("e", newgraphs,
                                                    ecount(res),
                                                    maps=maps)
    }
  } else {

    if (!keep.all.vertices) {
      minsize <- min(sapply(graphs, vcount))
      graphs <- lapply(graphs, function(g) {
        vc <- vcount(g)
        if (vc > minsize) {
          g <- g - (minsize+1):vc
        }
        g
      })
    }

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

#' Union of two or more sets
#'
#' This is an S3 generic function. See \code{methods("union")}
#' for the actual implementations for various S3 classes. Initially
#' it is implemented for igraph graphs and igraph vertex and edge
#' sequences. See
#' \code{\link{union.igraph}}, and
#' \code{\link{union.igraph.vs}}.
#'
#' @param ... Arguments, their number and interpretation depends on
#' the function that implements \code{union}.
#' @return Depends on the function that implements this method.
#'
#' @export

union <- function(...)
  UseMethod("union")

#' @method union default
#' @export

union.default <- function(...) {
  base::union(...)
}

#' Union of graphs
#' 
#' The union of two or more graphs are created. The graphs may have identical
#' or overlapping vertex sets.
#' 
#' \code{union} creates the union of two or more graphs.  Edges which are
#' included in at least one graph will be part of the new graph. This function
#' can be also used via the \%u\% operator.
#' 
#' If the \code{byname} argument is \code{TRUE} (or \code{auto} and all graphs
#' are named), then the operation is performed on symbolic vertex names instead
#' of the internal numeric vertex ids.
#' 
#' \code{union} keeps the attributes of all graphs. All graph, vertex and
#' edge attributes are copied to the result. If an attribute is present in
#' multiple graphs and would result a name clash, then this attribute is
#' renamed by adding suffixes: _1, _2, etc.
#' 
#' The \code{name} vertex attribute is treated specially if the operation is
#' performed based on symbolic vertex names. In this case \code{name} must be
#' present in all graphs, and it is not renamed in the result graph.
#' 
#' An error is generated if some input graphs are directed and others are
#' undirected.
#' 
#' @aliases graph.union %u%
#' @param \dots Graph objects or lists of graph objects.
#' @param byname A logical scalar, or the character scalar \code{auto}. Whether
#' to perform the operation based on symbolic vertex names. If it is
#' \code{auto}, that means \code{TRUE} if all graphs are named and \code{FALSE}
#' otherwise. A warning is generated if \code{auto} and some (but not all)
#' graphs are named.
#' @return A new graph object.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @method union igraph
#' @export
#' @keywords graphs
#' @examples
#' 
#' ## Union of two social networks with overlapping sets of actors
#' net1 <- graph_from_literal(D-A:B:F:G, A-C-F-A, B-E-G-B, A-B, F-G,
#'                   H-F:G, H-I-J)
#' net2 <- graph_from_literal(D-A:F:Y, B-A-X-F-H-Z, F-Y)
#' str(net1 %u% net2)

union.igraph <- function(..., byname="auto") {
  .igraph.graph.union.or.intersection("R_igraph_union", ..., byname=byname,
                                      keep.all.vertices=TRUE)
}

#' @export

"%u%" <- function(x,y) {
  union(x,y)
}

#' Intersection of two or more sets
#'
#' This is an S3 generic function. See \code{methods("intersection")}
#' for the actual implementations for various S3 classes. Initially
#' it is implemented for igraph graphs and igraph vertex and edge
#' sequences. See
#' \code{\link{intersection.igraph}}, and
#' \code{\link{intersection.igraph.vs}}.
#'
#' @param ... Arguments, their number and interpretation depends on
#' the function that implements \code{intersection}.
#' @return Depends on the function that implements this method.
#'
#' @export

intersection <- function(...)
  UseMethod("intersection")

#' Intersection of graphs
#' 
#' The intersection of two or more graphs are created.  The graphs may have
#' identical or overlapping vertex sets.
#' 
#' \code{intersection} creates the intersection of two or more graphs:
#' only edges present in all graphs will be included.  The corresponding
#' operator is \%s\%.
#' 
#' If the \code{byname} argument is \code{TRUE} (or \code{auto} and all graphs
#' are named), then the operation is performed on symbolic vertex names instead
#' of the internal numeric vertex ids.
#' 
#' \code{intersection} keeps the attributes of all graphs. All graph,
#' vertex and edge attributes are copied to the result. If an attribute is
#' present in multiple graphs and would result a name clash, then this
#' attribute is renamed by adding suffixes: _1, _2, etc.
#' 
#' The \code{name} vertex attribute is treated specially if the operation is
#' performed based on symbolic vertex names. In this case \code{name} must be
#' present in all graphs, and it is not renamed in the result graph.
#' 
#' An error is generated if some input graphs are directed and others are
#' undirected.
#' 
#' @aliases graph.intersection %s%
#' @param \dots Graph objects or lists of graph objects.
#' @param byname A logical scalar, or the character scalar \code{auto}. Whether
#' to perform the operation based on symbolic vertex names. If it is
#' \code{auto}, that means \code{TRUE} if all graphs are named and \code{FALSE}
#' otherwise. A warning is generated if \code{auto} and some (but not all)
#' graphs are named.
#' @param keep.all.vertices Logical scalar, whether to keep vertices that only
#' appear in a subset of the input graphs.
#' @return A new graph object.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @method intersection igraph
#' @export
#' @keywords graphs
#' @examples
#' 
#' ## Common part of two social networks
#' net1 <- graph_from_literal(D-A:B:F:G, A-C-F-A, B-E-G-B, A-B, F-G,
#'                   H-F:G, H-I-J)
#' net2 <- graph_from_literal(D-A:F:Y, B-A-X-F-H-Z, F-Y)
#' str(net1 %s% net2)

intersection.igraph <- function(..., byname="auto",
                                keep.all.vertices=TRUE) {
  .igraph.graph.union.or.intersection("R_igraph_intersection", ...,
                                      byname=byname,
                                      keep.all.vertices=keep.all.vertices)
}

#' @export

"%s%" <- function(x,y) {
  intersection(x,y)
}

#' Difference of two sets
#'
#' This is an S3 generic function. See \code{methods("difference")}
#' for the actual implementations for various S3 classes. Initially
#' it is implemented for igraph graphs (difference of edges in two graphs),
#' and igraph vertex and edge sequences. See
#' \code{\link{difference.igraph}}, and
#' \code{\link{difference.igraph.vs}}.
#'
#' @param ... Arguments, their number and interpretation depends on
#' the function that implements \code{difference}.
#' @return Depends on the function that implements this method.
#'
#' @export

difference <- function(...)
  UseMethod("difference")


#' Difference of graphs
#' 
#' The difference of two graphs are created.
#' 
#' \code{difference} creates the difference of two graphs. Only edges
#' present in the first graph but not in the second will be be included in the
#' new graph. The corresponding operator is \%m\%.
#' 
#' If the \code{byname} argument is \code{TRUE} (or \code{auto} and the graphs
#' are all named), then the operation is performed based on symbolic vertex
#' names. Otherwise numeric vertex ids are used.
#' 
#' \code{difference} keeps all attributes (graph, vertex and edge) of the
#' first graph.
#' 
#' Note that \code{big} and \code{small} must both be directed or both be
#' undirected, otherwise an error message is given.
#' 
#' @aliases graph.difference %m%
#' @param big The left hand side argument of the minus operator. A directed or
#' undirected graph.
#' @param small The right hand side argument of the minus operator. A directed
#' ot undirected graph.
#' @param byname A logical scalar, or the character scalar \code{auto}. Whether
#' to perform the operation based on symbolic vertex names. If it is
#' \code{auto}, that means \code{TRUE} if both graphs are named and
#' \code{FALSE} otherwise. A warning is generated if \code{auto} and one graph,
#' but not both graphs are named.
#' @param ... Ignored, included for S3 compatibility.
#' @return A new graph object.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @method difference igraph
#' @export
#' @keywords graphs
#' @examples
#' 
#' ## Create a wheel graph
#' wheel <- union(make_ring(10),
#'                      make_star(11, center=11, mode="undirected"))
#' V(wheel)$name <- letters[seq_len(vcount(wheel))]
#' 
#' ## Subtract a star graph from it
#' sstar <- make_star(6, center=6, mode="undirected")
#' V(sstar)$name <- letters[c(1,3,5,7,9,11)]
#' G <- wheel %m% sstar
#' str(G)
#' plot(G, layout=layout_nicely(wheel))

difference.igraph <- function(big, small, byname="auto", ...) {

  if (!is_igraph(big) || !is_igraph(small)) {
    stop("argument is not a graph")
  }
  if (byname != "auto" && !is.logical(byname)) {
    stop("`bynam' must be \"auto\", or logical")
  }
  nonamed <- is_named(big) + is_named(small)
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
    big <- permute(big, perm)

    on.exit(.Call("R_igraph_finalizer", PACKAGE="igraph"))
    res <- .Call("R_igraph_difference", big, small,
                 PACKAGE="igraph")
    permute(res, match(V(res)$name, bnames))

  } else {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    .Call("R_igraph_difference", big, small,
          PACKAGE="igraph")
  }
}

#' @export

"%m%" <- function(x,y) {
  difference(x,y)
}



#' Complementer of a graph
#' 
#' A complementer graph contains all edges that were not present in the input
#' graph.
#' 
#' \code{complementer} creates the complementer of a graph. Only edges
#' which are \emph{not} present in the original graph will be included in the
#' new graph.
#' 
#' \code{complementer} keeps graph and vertex attriubutes, edge
#' attributes are lost.
#'
#' @aliases graph.complementer
#' @param graph The input graph, can be directed or undirected.
#' @param loops Logical constant, whether to generate loop edges.
#' @return A new graph object.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @export
#' @keywords graphs
#' @examples
#' 
#' ## Complementer of a ring
#' g <- make_ring(10)
#' complementer(g)
#' 
#' ## A graph and its complementer give together the full graph
#' g <- make_ring(10)
#' gc <- complementer(g)
#' gu <- union(g, gc)
#' gu
#' graph.isomorphic(gu, make_full_graph(vcount(g)))
#' 
complementer <- function(graph, loops=FALSE) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_complementer", graph, as.logical(loops),
        PACKAGE="igraph")
}



#' Compose two graphs as binary relations
#' 
#' Relational composition of two graph.
#' 
#' \code{compose} creates the relational composition of two graphs. The
#' new graph will contain an (a,b) edge only if there is a vertex c, such that
#' edge (a,c) is included in the first graph and (c,b) is included in the
#' second graph. The corresponding operator is \%c\%.
#' 
#' The function gives an error if one of the input graphs is directed and the
#' other is undirected.
#' 
#' If the \code{byname} argument is \code{TRUE} (or \code{auto} and the graphs
#' are all named), then the operation is performed based on symbolic vertex
#' names. Otherwise numeric vertex ids are used.
#' 
#' \code{compose} keeps the attributes of both graphs. All graph, vertex
#' and edge attributes are copied to the result. If an attribute is present in
#' multiple graphs and would result a name clash, then this attribute is
#' renamed by adding suffixes: _1, _2, etc.
#' 
#' The \code{name} vertex attribute is treated specially if the operation is
#' performed based on symbolic vertex names. In this case \code{name} must be
#' present in both graphs, and it is not renamed in the result graph.
#' 
#' Note that an edge in the result graph corresponds to two edges in the input,
#' one in the first graph, one in the second. This mapping is not injective and
#' several edges in the result might correspond to the same edge in the first
#' (and/or the second) graph. The edge attributes in the result graph are
#' updated accordingly.
#' 
#' Also note that the function may generate multigraphs, if there are more than
#' one way to find edges (a,b) in g1 and (b,c) in g2 for an edge (a,c) in the
#' result. See \code{\link{simplify}} if you want to get rid of the multiple
#' edges.
#' 
#' The function may create loop edges, if edges (a,b) and (b,a) are present in
#' g1 and g2, respectively, then (a,a) is included in the result. See
#' \code{\link{simplify}} if you want to get rid of the self-loops.
#' 
#' @aliases graph.compose %c%
#' @param g1 The first input graph.
#' @param g2 The second input graph.
#' @param byname A logical scalar, or the character scalar \code{auto}. Whether
#' to perform the operation based on symbolic vertex names. If it is
#' \code{auto}, that means \code{TRUE} if both graphs are named and
#' \code{FALSE} otherwise. A warning is generated if \code{auto} and one graph,
#' but not both graphs are named.
#' @return A new graph object.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @export
#' @keywords graphs
#' @examples
#' 
#' g1 <- make_ring(10)
#' g2 <- make_star(10, mode="undirected")
#' gc <- compose(g1, g2)
#' str(gc)
#' str(simplify(gc))
#' 
compose <- function(g1, g2, byname="auto") {

  if (!is_igraph(g1) || !is_igraph(g2)) {
    stop("Not a graph object")
  }

  if (byname != "auto" && !is.logical(byname)) {
    stop("`byname' must be \"auto\", or logical")
  }
  nonamed <- is_named(g1) + is_named(g2)
  if (byname == "auto") {
    byname <- nonamed == 2
    if (nonamed == 1) {
      warning("One, but not both graphs are named, not using vertex names")
    }
  } else if (byname && nonamed != 2) {
    stop("Some graphs are not named")
  }

  if (byname) {
    uninames <- unique(c(V(g1)$name, V(g2)$name))
    if (vcount(g1) < length(uninames)) {
      g1 <- g1 + setdiff(uninames, V(g1)$name)
    }
    if (vcount(g2) < length(uninames)) {
      g2 <- g2 + setdiff(uninames, V(g2)$name)
    }
    if (any(uninames != V(g1)$name)) {
      g1 <- permute(g1, match(V(g1)$name, uninames))
    }
    if (any(uninames != V(g2)$name)) {
      g2 <- permute(g2, match(V(g2)$name, uninames))
    }
  }

  edgemaps <- (length(edge_attr_names(g1)) != 0 ||
               length(edge_attr_names(g2)) != 0)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_compose", g1, g2, edgemaps,
               PACKAGE="igraph")
  maps <- list(res$edge_map1, res$edge_map2)
  res <- res$graph

  ## We might need to rename all attributes
  graphs <- list(g1, g2)
  graph.attributes(res) <- rename.attr.if.needed("g", graphs)

  if (byname) {
    vertex.attributes(res) <-
      rename.attr.if.needed("v", graphs, vcount(res), ignore="name")
    V(res)$name <- uninames
  } else {
    vertex.attributes(res) <- rename.attr.if.needed("v", graphs,
                                                    vcount(res))
  }

  if (edgemaps) {
    edge.attributes(res) <- rename.attr.if.needed("e", graphs, ecount(res),
                                                  maps2=maps)
  }

  res
}

#' @export

"%c%" <- function(x,y) {
  compose(x,y)
}

#' @export

edge <- function(...) {
  structure(list(...), class="igraph.edge")
}

#' @export

edges <- edge

#' @export

vertex <- function(...) {
  structure(list(...), class="igraph.vertex")
}

#' @export

vertices <- vertex

#' @export

path <- function(...) {
  structure(list(...), class="igraph.path")
}

#' @method "+" igraph
#' @export

`+.igraph` <- function(e1, e2) {
  if (!is_igraph(e1) && is_igraph(e2)) {
    tmp <- e1
    e1 <- e2
    e2 <- tmp
  }
  if (is_igraph(e2) && is_named(e1) && is_named(e2)) {
    ## Union of graphs
    res <- union(e1, e2)
  } else if (is_igraph(e2)) {
    ## Disjoint union of graphs
    res <- disjoint_union(e1,e2)

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
    res <- add_edges(e1, as.igraph.vs(e1, toadd), attr=attr)

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
    res <- add_vertices(e1, la, attr=e2)

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
      res <- add_edges(e1, toadd, attr=attr)
    } else {
      res <- e1
    }
    
  } else if (is.numeric(e2) && length(e2)==1) {
    ## Adding some isolate vertices
    res <- add_vertices(e1, e2)

  } else if (is.character(e2)) {
    ## Adding named vertices
    res <- add_vertices(e1, length(e2), name=e2)
    
  } else {
    stop("Cannot add unknown type to igraph graph")
  }
  res
}

#' @method "-" igraph
#' @export
 
`-.igraph` <- function(e1, e2) {
  if (missing(e2)) {
    stop("Non-numeric argument to negation operator")
  }
  if (is_igraph(e2)) {
    res <- difference(e1, e2)
  } else if ("igraph.vertex" %in% class(e2)) {
    res <- delete_vertices(e1, unlist(e2, recursive=FALSE))
  } else if ("igraph.edge" %in% class(e2)) {
    res <- delete_edges(e1, unlist(e2, recursive=FALSE))
  } else if ("igraph.path" %in% class(e2)) {
    todel <- unlist(e2, recursive=FALSE)
    lt <- length(todel)
    if (lt >= 2) {
      todel <- paste(todel[-lt], todel[-1], sep="|")
      res <- delete_edges(e1, todel)
    } else {
      res <- e1
    }
  } else if ("igraph.vs" %in% class(e2)) {
    res <- delete_vertices(e1, e2)
  } else if ("igraph.es" %in% class(e2)) {
    res <- delete_edges(e1, e2)
  } else if (is.numeric(e2) || is.character(e2)) {
    res <- delete_vertices(e1, e2)
  } else {
    stop("Cannot substract unknown type from igraph graph")
  }
  res
}
