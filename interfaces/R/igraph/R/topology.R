
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

graph.get.isomorphisms.vf2 <- function(graph1, graph2, vertex.color1,
                                       vertex.color2, edge.color1,
                                       edge.color2) {
  # Argument checks
  if (!is.igraph(graph1)) { stop("Not a graph object") }
  if (!is.igraph(graph2)) { stop("Not a graph object") }
  if (missing(vertex.color1)) { 
    if ("color" %in% list.vertex.attributes(graph1)) { 
      vertex.color1 <- V(graph1)$color 
    } else { 
      vertex.color1 <- NULL 
    } 
  } 
  if (!is.null(vertex.color1)) { 
    vertex.color1 <- as.integer(vertex.color1)-1L 
  }
  if (missing(vertex.color2)) { 
    if ("color" %in% list.vertex.attributes(graph2)) { 
      vertex.color2 <- V(graph2)$color 
    } else { 
      vertex.color2 <- NULL 
    } 
  } 
  if (!is.null(vertex.color2)) { 
    vertex.color2 <- as.integer(vertex.color2)-1L 
  }
  if (missing(edge.color1)) { 
    if ("color" %in% list.edge.attributes(graph1)) { 
      edge.color1 <- E(graph1)$color 
    } else { 
      edge.color1 <- NULL 
    } 
  } 
  if (!is.null(edge.color1)) { 
    edge.color1 <- as.integer(edge.color1)-1L 
  }
  if (missing(edge.color2)) { 
    if ("color" %in% list.edge.attributes(graph2)) { 
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
  
  res <- lapply(res, "+", 1)

  if (getIgraphOpt("add.vertex.names") && is.named(graph2)) {
    for (i in seq_along(res)) {
      names(res[[i]]) <- V(graph2)$name[ res[[i]] ]
    }
  }
  res
}

graph.get.subisomorphisms.vf2 <- function(graph1, graph2, vertex.color1,
                                          vertex.color2, edge.color1,
                                          edge.color2) {
  # Argument checks
  if (!is.igraph(graph1)) { stop("Not a graph object") }
  if (!is.igraph(graph2)) { stop("Not a graph object") }
  if (missing(vertex.color1)) { 
    if ("color" %in% list.vertex.attributes(graph1)) { 
      vertex.color1 <- V(graph1)$color 
    } else { 
      vertex.color1 <- NULL 
    } 
  } 
  if (!is.null(vertex.color1)) { 
    vertex.color1 <- as.integer(vertex.color1)-1L 
  }
  if (missing(vertex.color2)) { 
    if ("color" %in% list.vertex.attributes(graph2)) { 
      vertex.color2 <- V(graph2)$color 
    } else { 
      vertex.color2 <- NULL 
    } 
  } 
  if (!is.null(vertex.color2)) { 
    vertex.color2 <- as.integer(vertex.color2)-1L 
  }
  if (missing(edge.color1)) { 
    if ("color" %in% list.edge.attributes(graph1)) { 
      edge.color1 <- E(graph1)$color 
    } else { 
      edge.color1 <- NULL 
    } 
  } 
  if (!is.null(edge.color1)) { 
    edge.color1 <- as.integer(edge.color1)-1L 
  }
  if (missing(edge.color2)) { 
    if ("color" %in% list.edge.attributes(graph2)) { 
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

  res <- lapply(res, "+", 1)

  if (getIgraphOpt("add.vertex.names") && is.named(graph2)) {
    for (i in seq_along(res)) {
      names(res[[i]]) <- V(graph2)$name[ res[[i]] ]
    }
  }
  res
}

graph.isoclass.subgraph <- function(graph, vids) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  vids <- as.igraph.vs(graph, vids)-1

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_isoclass_subgraph", graph, vids,
        PACKAGE="igraph")
  res
}

graph.subisomorphic.lad <- function(pattern, target, domains=NULL,
                                    induced=FALSE, map=TRUE, all.maps=FALSE,
                                    time.limit=Inf) {
  # Argument checks
  if (!is.igraph(pattern)) { stop("Not a graph object") }
  if (!is.igraph(target)) { stop("Not a graph object") }
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
    if (getIgraphOpt("add.vertex.names") && is.named(target)) {
      names(res$map) <- V(target)$name[res$map]
    }
  }
  if (all.maps) {
    res$maps <- lapply(res$maps, function(x) x + 1)
    if (getIgraphOpt("add.vertex.names") && is.named(target)) {
      for (i in seq_along(res$maps)) {
        names(res$maps[[i]]) <- V(target)$name[ res$maps[[i]] ]
      }
    }
  }
  res
}


# Rest generated by stimulus

