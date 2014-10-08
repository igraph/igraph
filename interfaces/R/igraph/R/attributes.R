#   IGraph R package
#   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

##
## The brand new attribute interface:
##
## g(graph)$name                   # get a graph attribute
## g(graph)$name <- "Ring"         # set a graph attribute
## 
## v(graph)$color <- "red"         # set vertex attribute
## v(graph)$color[1:5] <- "blue"   
## v(graph)$color[c(5,6,7)]        # get vertex attribute
##
## e(graph)$weight <- 1            # set edge attribute
## e(graph)$weight[1:10]           # get edge attribute
##

#' @export

graph_attr <- function(graph, name) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  if (missing(name)) {
    graph.attributes(graph)
  } else {
    .Call("R_igraph_mybracket2", graph, 9L, 2L,
          PACKAGE="igraph")[[as.character(name)]]
  }
}

#' @export

`graph_attr<-` <- function(graph, name, value) {
  if (missing(name)) {
    `graph.attributes<-`(graph, value)
  } else {
    set_graph_attr(graph, name, value)
  }
}

#' @export

set_graph_attr <- function(graph, name, value) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  .Call("R_igraph_mybracket3_set", graph, 9L, 2L, name, value,
        PACKAGE="igraph")
}

#' @export

graph.attributes <- function(graph) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_mybracket2_copy", graph, 9L, 2L,
        PACKAGE="igraph")
} 

#' @export

"graph.attributes<-" <- function(graph, value) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  if (!is.list(value) || (length(value) > 0 && is.null(names(value))) ||
      any(names(value) == "") || any(duplicated(names(value)))) {
    stop("Value must be a named list with unique names")
  }
            
  .Call("R_igraph_mybracket2_set", graph, 9L, 2L, value,
        PACKAGE="igraph")
}

#' @export

vertex_attr <- function(graph, name, index=V(graph)) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  if (missing(name)) {
    vertex.attributes(graph, index = index)
  } else {
    index <- as.igraph.vs(graph, index)
    myattr <- .Call("R_igraph_mybracket2", graph, 9L, 3L,
                    PACKAGE="igraph")[[as.character(name)]]
    myattr[index]
  }
}

#' @export

`vertex_attr<-` <- function(graph, name, index = V(graph), value) {
  if (missing(name)) {
    `vertex.attributes<-`(graph, index = index, value = value)
  } else {
    set_vertex_attr(graph, name = name, index = index, value = value)
  }
}

#' @export

set_vertex_attr <- function(graph, name, index=V(graph), value) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  single <- "single" %in% names(attributes(index)) && attr(index, "single")
  if (!missing(index)) { index <- as.igraph.vs(graph, index) }
  name <- as.character(name)
  vc <- vcount(graph)

  vattrs <- .Call("R_igraph_mybracket2", graph, 9L, 3L, PACKAGE="igraph")
  if (single) {
    vattrs[[name]][[index]] <- value
  } else {
    vattrs[[name]][index] <- value
  }
  length(vattrs[[name]]) <- vc
  
  .Call("R_igraph_mybracket2_set", graph, 9L, 3L, vattrs, PACKAGE="igraph")
}

#' @export

vertex.attributes <- function(graph, index = V(graph)) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  if (!missing(index)) {
    index <- as.igraph.vs(graph, index)
  }

  res <- .Call("R_igraph_mybracket2_copy", graph, 9L, 3L, PACKAGE="igraph")

  if (!missing(index) &&
      (length(index) != vcount(graph) || any(index != V(graph)))) {
    for (i in seq_along(value)) {
      value[[i]] <- value[[i]][index]
    }
  }
  res
}

#' @export

"vertex.attributes<-" <- function(graph, index = V(graph), value) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  if (!is.list(value) || (length(value) > 0 && is.null(names(value))) ||
      any(names(value) == "") || any(duplicated(names(value)))) {
    stop("Value must be a named list with unique names")
  }
  if ( any(sapply(value, length) != length(index)) ) {
    stop("Invalid attribute value length, must match number of vertices")
  }

  if (!missing(index)) {
    index <- as.igraph.vs(graph, index)
    if (any(duplicated(index)) || anyNA(index)) {
      stop("Invalid vertices in index")
    }
  }

  if (!missing(index) &&
      (length(index) != vcount(graph) || any(index != V(graph)))) {
    vs <- V(graph)
    for (i in seq_along(value)) {
      tmp <- value[[i]]
      length(tmp) <- 0
      length(tmp) <- length(vs)
      tmp[index] <- value[[i]]
      value[[i]] <- tmp
    }
  }

  .Call("R_igraph_mybracket2_set", graph, 9L, 3L, value,
        PACKAGE="igraph")
}

#' @export

edge_attr <- function(graph, name, index=E(graph)) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  if (missing(name)) {
    edge.attributes(graph, name)
  } else {
    name <- as.character(name)
    index <- as.igraph.es(graph, index)
    myattr <- .Call("R_igraph_mybracket2", graph, 9L, 4L,
                    PACKAGE="igraph")[[name]]
    myattr[index]
  }
}

#' @export

`edge_attr<-` <- function(graph, name, index = E(graph), value) {
  if (missing(name)) {
    `edge.attributes<-`(graph, index = index, value = value)
  } else {
    set_edge_attr(graph, name = name, index = index, value = value)
  }
}

#' @export

set_edge_attr <- function(graph, name, index=E(graph), value) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  single <- "single" %in% names(attributes(index)) && attr(index, "single")
  name <- as.character(name)
  index <- as.igraph.es(graph, index)
  ec <- ecount(graph)

  eattrs <- .Call("R_igraph_mybracket2", graph, 9L, 4L, PACKAGE="igraph")
  if (single) {
    eattrs[[name]][[index]] <- value
  } else {
    eattrs[[name]][index] <- value
  }
  length(eattrs[[name]]) <- ec

  .Call("R_igraph_mybracket2_set", graph, 9L, 4L, eattrs, PACKAGE="igraph")
}

#' @export

edge.attributes <- function(graph, index = E(graph)) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  if (!missing(index)) {
    index <- as.igraph.es(graph, index)
  }

  res <- .Call("R_igraph_mybracket2_copy", graph, 9L, 4L, PACKAGE="igraph")

  if (!missing(index) &&
      (length(index) != ecount(graph) || any(index != E(graph)))) {
    for (i in seq_along(value)) {
      value[[i]] <- value[[i]][index]
    }
  }
  res
}

#' @export

"edge.attributes<-" <- function(graph, index = E(graph), value) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  if (!is.list(value) || (length(value) > 0 && is.null(names(value))) ||
      any(names(value) == "") || any(duplicated(names(value)))) {
    stop("Value must be a named list with unique names")
  }
  if ( any(sapply(value, length) != length(index)) ) {
    stop("Invalid attribute value length, must match number of edges")
  }

  if (!missing(index)) {
    index <- as.igraph.es(graph, index)
    if (any(duplicated(index)) || anyNA(index)) {
      stop("Invalid edges in index")
    }
  }

  if (!missing(index) &&
      (length(index) != ecount(graph) || any(index != E(graph)))) {
    es <- E(graph)
    for (i in seq_along(value)) {
      tmp <- value[[i]]
      length(tmp) <- 0
      length(tmp) <- length(es)
      tmp[index] <- value[[i]]
      value[[i]] <- tmp
    }
  }
  
  .Call("R_igraph_mybracket2_set", graph, 9L, 4L, value,
        PACKAGE="igraph")
}

#' @export

graph_attr_names <- function(graph) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  res <- .Call("R_igraph_mybracket2_names", graph, 9L, 2L, PACKAGE="igraph")
  if (is.null(res)) { res <- character() }
  res
}

#' @export

vertex_attr_names <- function(graph) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  res <- .Call("R_igraph_mybracket2_names", graph, 9L, 3L, PACKAGE="igraph")
                     
  if (is.null(res)) { res <- character() }
  res
}

#' @export

edge_attr_names <- function(graph) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  res <- .Call("R_igraph_mybracket2_names", graph, 9L, 4L, PACKAGE="igraph")
  if (is.null(res)) { res <- character() }
  res
}

#' @export

delete_graph_attr <- function(graph, name) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  name <- as.character(name)
  if (!name %in% graph_attr_names(graph)) {
    stop("No such graph attribute: ", name)
  }

  gattr <- .Call("R_igraph_mybracket2", graph, 9L, 2L, PACKAGE="igraph")
  gattr[[name]] <- NULL
  
  .Call("R_igraph_mybracket2_set", graph, 9L, 2L, gattr, PACKAGE="igraph")
}

#' @export

delete_vertex_attr <- function(graph, name) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  name <- as.character(name)
  if (!name %in% vertex_attr_names(graph)) {
    stop("No such vertex attribute: ", name)
  }

  vattr <- .Call("R_igraph_mybracket2", graph, 9L, 3L, PACKAGE="igraph")
  vattr[[name]] <- NULL
  
  .Call("R_igraph_mybracket2_set", graph, 9L, 3L, vattr, PACKAGE="igraph")
}

#' @export

delete_edge_attr <- function(graph, name) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  name <- as.character(name)
  if (!name %in% edge_attr_names(graph)) {
    stop("No such edge attribute: ", name)
  }

  eattr <- .Call("R_igraph_mybracket2", graph, 9L, 4L, PACKAGE="igraph")
  eattr[[name]] <- NULL
  
  .Call("R_igraph_mybracket2_set", graph, 9L, 4L, eattr, PACKAGE="igraph")
}

#############



#' Named graphs
#' 
#' An igraph graph is named, if there is a symbolic name associated with its
#' vertices.
#' 
#' In igraph vertices can always be identified and specified via their numeric
#' vertex ids. This is, however, not always convenient, and in many cases there
#' exist symbolic ids that correspond to the vertices. To allow this more
#' flexible identification of vertices, one can assign a vertex attribute
#' called \sQuote{name} to an igraph graph. After doing this, the symbolic
#' vertex names can be used in all igraph functions, instead of the numeric
#' ids.
#' 
#' Note that the uniqueness of vertex names are currently not enforced in
#' igraph, you have to check that for yourself, when assigning the vertex
#' names.
#'
#' @aliases is.named
#' @param graph The input graph.
#' @return A logical scalar.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- ring(10)
#' is_named(g)
#' V(g)$name <- letters[1:10]
#' is_named(g)
#' neighbors(g, "a")
#' 
is_named <- function(graph) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  "name" %in% vertex_attr_names(graph)
}



#' Weighted graphs
#' 
#' In weighted graphs, a real number is assigned to each (directed or
#' undirected) edge.
#' 
#' In igraph edge weights are represented via an edge attribute, called
#' \sQuote{weight}. The \code{is_weighted} function only checks that such an
#' attribute exists. (It does not even checks that it is a numeric edge
#' attribute.)
#' 
#' Edge weights are used for different purposes by the different functions.
#' E.g. shortest path functions use it as the cost of the path; community
#' finding methods use it as the strength of the relationship between two
#' vertices, etc. Check the manual pages of the functions working with weighted
#' graphs for details.
#'
#' @aliases is.weighted
#' @param graph The input graph.
#' @return A logical scalar.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- ring(10)
#' get.shortest.paths(g, 8, 2)
#' E(g)$weight <- seq_len(ecount(g))
#' get.shortest.paths(g, 8, 2)
#' 
is_weighted <- function(graph) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  "weight" %in% edge_attr_names(graph)
}

#' @rdname bipartite_graph
#' @export

is_bipartite <- function(graph) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  "type" %in% vertex_attr_names(graph)
}

#############

igraph.i.attribute.combination <- function(comb) {
  if (is.function(comb)) {
    comb <- list(comb)
  }
  comb <- as.list(comb)
  if (any(!sapply(comb, function(x)
                  is.function(x) || (is.character(x) && length(x)==1)))) {
    stop("Attribute combination element must be a function or character scalar")
  }
  if (is.null(names(comb))) {
    names(comb) <- rep("", length(comb))
  }
  if (any(duplicated(names(comb)))) {
    warning("Some attributes are duplicated")
  }
  comb <- lapply(comb, function(x) {
    if (!is.character(x)) {
      x
    } else {
      known <- data.frame(n=c("ignore", "sum", "prod", "min", "max", "random",
                            "first", "last", "mean", "median", "concat"),
                          i=c(0,3,4,5,6,7,8,9,10,11,12), stringsAsFactors=FALSE)
      x <- pmatch(tolower(x), known[,1])
      if (is.na(x)) {
        stop("Unknown/unambigous attribute combination specification")
      }
      known[,2][x]
    }
  })

  comb
}
