
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

#' @export

graph.get.isomorphisms.vf2 <- function(graph1, graph2, vertex.color1,
                                       vertex.color2, edge.color1,
                                       edge.color2) {
  # Argument checks
  if (!is_igraph(graph1)) { stop("Not a graph object") }
  if (!is_igraph(graph2)) { stop("Not a graph object") }
  if (missing(vertex.color1)) { 
    if ("color" %in% vertex_attr_names(graph1)) { 
      vertex.color1 <- V(graph1)$color 
    } else { 
      vertex.color1 <- NULL 
    } 
  } 
  if (!is.null(vertex.color1)) { 
    vertex.color1 <- as.integer(vertex.color1)-1L 
  }
  if (missing(vertex.color2)) { 
    if ("color" %in% vertex_attr_names(graph2)) { 
      vertex.color2 <- V(graph2)$color 
    } else { 
      vertex.color2 <- NULL 
    } 
  } 
  if (!is.null(vertex.color2)) { 
    vertex.color2 <- as.integer(vertex.color2)-1L 
  }
  if (missing(edge.color1)) { 
    if ("color" %in% edge_attr_names(graph1)) { 
      edge.color1 <- E(graph1)$color 
    } else { 
      edge.color1 <- NULL 
    } 
  } 
  if (!is.null(edge.color1)) { 
    edge.color1 <- as.integer(edge.color1)-1L 
  }
  if (missing(edge.color2)) { 
    if ("color" %in% edge_attr_names(graph2)) { 
      edge.color2 <- E(graph2)$color 
    } else { 
      edge.color2 <- NULL 
    } 
  } 
  if (!is.null(edge.color2)) { 
    edge.color2 <- as.integer(edge.color2)-1L 
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_get_isomorphisms_vf2", graph1, graph2, vertex.color1,
               vertex.color2, edge.color1, edge.color2,
               PACKAGE="igraph")
  
  lapply(res, function(x) V(graph2)[x + 1])
}

#' @export

graph.get.subisomorphisms.vf2 <- function(graph1, graph2, vertex.color1,
                                          vertex.color2, edge.color1,
                                          edge.color2) {
  # Argument checks
  if (!is_igraph(graph1)) { stop("Not a graph object") }
  if (!is_igraph(graph2)) { stop("Not a graph object") }
  if (missing(vertex.color1)) { 
    if ("color" %in% vertex_attr_names(graph1)) { 
      vertex.color1 <- V(graph1)$color 
    } else { 
      vertex.color1 <- NULL 
    } 
  } 
  if (!is.null(vertex.color1)) { 
    vertex.color1 <- as.integer(vertex.color1)-1L 
  }
  if (missing(vertex.color2)) { 
    if ("color" %in% vertex_attr_names(graph2)) { 
      vertex.color2 <- V(graph2)$color 
    } else { 
      vertex.color2 <- NULL 
    } 
  } 
  if (!is.null(vertex.color2)) { 
    vertex.color2 <- as.integer(vertex.color2)-1L 
  }
  if (missing(edge.color1)) { 
    if ("color" %in% edge_attr_names(graph1)) { 
      edge.color1 <- E(graph1)$color 
    } else { 
      edge.color1 <- NULL 
    } 
  } 
  if (!is.null(edge.color1)) { 
    edge.color1 <- as.integer(edge.color1)-1L 
  }
  if (missing(edge.color2)) { 
    if ("color" %in% edge_attr_names(graph2)) { 
      edge.color2 <- E(graph2)$color 
    } else { 
      edge.color2 <- NULL 
    } 
  } 
  if (!is.null(edge.color2)) { 
    edge.color2 <- as.integer(edge.color2)-1L 
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_get_subisomorphisms_vf2", graph1, graph2,
               vertex.color1, vertex.color2, edge.color1, edge.color2,
               PACKAGE="igraph")

  lapply(res, function(x) V(graph1)[x + 1])
}

#' @export

graph.isoclass.subgraph <- function(graph, vids) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  vids <- as.igraph.vs(graph, vids)-1

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_isoclass_subgraph", graph, vids,
        PACKAGE="igraph")
  res
}

#' @export

graph.subisomorphic.lad <- function(pattern, target, domains=NULL,
                                    induced=FALSE, map=TRUE, all.maps=FALSE,
                                    time.limit=Inf) {
  # Argument checks
  if (!is_igraph(pattern)) { stop("Not a graph object") }
  if (!is_igraph(target)) { stop("Not a graph object") }
  induced <- as.logical(induced)
  if (time.limit==Inf) {
    time.limit <- 0L
  } else {
    time.limit <- as.integer(time.limit)
  }
  map <- as.logical(map)
  all.maps <- as.logical(all.maps)
  if (!is.null(domains)) {
    if (!is.list(domains)) {
      stop("`domains' must be a list of vertex vectors from `target'")
    }
    if (length(domains) != vcount(pattern)) {
      stop("`domains' length and `pattern' number of vertices must match")
    }
    domains <- lapply(domains, function(x) as.igraph.vs(target, x)-1)
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_subisomorphic_lad", pattern, target, domains,
               induced, time.limit, map, all.maps,
               PACKAGE="igraph")

  if (map) {
    res$map <- res$map + 1
    if (igraph_opt("add.vertex.names") && is_named(target)) {
      names(res$map) <- V(target)$name[res$map]
    }
  }
  if (all.maps) res$maps <- lapply(res$maps, function(x) V(target)[x+1])

  res
}

## ----------------------------------------------------------------------
## NEW API

#' @export

isomorphic <- function(graph1, graph2, method = c("auto", "direct",
                 "vf2", "bliss"), ...) {

  if (!is_igraph(graph1)) { stop("Not a graph object") }
  if (!is_igraph(graph2)) { stop("Not a graph object") }
  method <- igraph.match.arg(method)

  if (method == "auto") {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    .Call("R_igraph_isomorphic", graph1, graph2, PACKAGE="igraph")
    
  } else if (method == "direct") {
    on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
    .Call("R_igraph_isomorphic_34", graph1, graph2, PACKAGE = "igraph")

  } else if (method == "vf2") {
    graph.isomorphic.vf2(graph1, graph2, ...)$iso
    
  } else if (method == "bliss") {
    graph.isomorphic.bliss(graph1, graph2, ...)$iso
    
  }
}

#' @export

is_isomorphic_to <- isomorphic

#' @export

subgraph_isomorphic <- function(pattern, target,
                                method = c("auto", "lad", "vf2"), ...) {

  method <- igraph.match.arg(method)

  if (method == "auto") method <- "lad"

  if (method == "lad") {
    graph.subisomorphic.lad(pattern, target, map = FALSE, all.maps = FALSE,
                            ...)$iso

  } else if (method == "vf2") {
    graph.subisomorphic.vf2(target, pattern, ...)$iso
    
  }
}


#' @export

is_subgraph_isomorphic_to <- subgraph_isomorphic


#' @include auto.R
#' @export

count_isomorphisms <- function(graph1, graph2, method = "vf2", ...) {

  method <- igraph.match.arg(method)

  if (method == "vf2") {
    graph.count.isomorphisms.vf2(graph1, graph2, ...)
  }

}


#' @export

count_subgraph_isomorphisms <- function(pattern, target,
                                        method = c("lad", "vf2"), ...) {

  method <- igraph.match.arg(method)

  if (method == "lad") {
    length(graph.subisomorphic.lad(pattern, target, all.maps = TRUE, ...)$maps)
    
  } else if (method == "vf2") {
    graph.count.subisomorphisms.vf2(target, pattern, ...)
  }

}


#' @export

isomorphisms <- function(graph1, graph2, method = "vf2", ...) {

  method <- igraph.match.arg(method)

  if (method == "vf2") {
    graph.get.isomorphisms.vf2(graph1, graph2, ...)
  }

}


#' @export

subgraph_isomorphisms <- function(pattern, target,
                                  method = c("lad", "vf2"), ...) {

  method <- igraph.match.arg(method)

  if (method == "lad") {
    graph.subisomorphic.lad(pattern, target, all.maps = TRUE, ...)$maps

  } else if (method == "vf2") {
    graph.get.subisomorphisms.vf2(target, pattern, ...)
  }

}


#' @export

isomorphism_class <- function(graph, v) {

  if (missing(v)) {
    graph.isoclass(graph)

  } else {
    graph.isoclass.subgraph(graph, v)
  }
  
}
