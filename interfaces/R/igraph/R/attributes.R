
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

get.graph.attribute <- function(graph, name) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_mybracket2", graph, 9L, 2L,
        PACKAGE="igraph")[[as.character(name)]]
}

set.graph.attribute <- function(graph, name, value) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  ## Trick to make R copy the graph
  newgraph <- graph
  attr(newgraph, "foo") <- NULL
  
  ## !!! Modifies the graph is place
  .Call("R_igraph_mybracket3_set", newgraph, 9L, 2L, name, value,
        PACKAGE="igraph")

  newgraph
}

graph.attributes <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_mybracket2_copy", graph, 9L, 2L,
        PACKAGE="igraph")
} 

"graph.attributes<-" <- function(graph, value) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (!is.list(value) || (length(value) > 0 && is.null(names(value))) ||
      any(names(value) == "") || any(duplicated(names(value)))) {
    stop("Value must be a named list with unique names")
  }
            
  ## Trick to make R copy the graph
  newgraph <- graph
  attr(newgraph, "foo") <- NULL

  ## !!! Modifies the graph in place
  .Call("R_igraph_mybracket2_set", newgraph, 9L, 2L, value,
        PACKAGE="igraph")
  
  newgraph
}

get.vertex.attribute <- function(graph, name, index=V(graph)) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  index <- as.igraph.vs(graph, index)  
  myattr <- .Call("R_igraph_mybracket2", graph, 9L, 3L,
                  PACKAGE="igraph")[[as.character(name)]]
  if (is.list(myattr) && length(index)==1) {
    myattr[[index]]
  } else {
    myattr[index]
  }
}

set.vertex.attribute <- function(graph, name, index=V(graph), value) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (!missing(index)) { index <- as.igraph.vs(graph, index) }
  name <- as.character(name)
  vc <- vcount(graph)

  vattrs <- .Call("R_igraph_mybracket2", graph, 9L, 3L, PACKAGE="igraph")
  vattrs[[name]][index] <- value
  length(vattrs[[name]]) <- vc
  
  ## Trick to make R copy the graph
  newgraph <- graph
  attr(newgraph, "foo") <- NULL

  ## !!! Modifies the graph in place
  .Call("R_igraph_mybracket2_set", newgraph, 9L, 3L, vattrs, PACKAGE="igraph")

  newgraph
}

vertex.attributes <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_mybracket2_copy", graph, 9L, 3L, PACKAGE="igraph")  
}

"vertex.attributes<-" <- function(graph, value) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (!is.list(value) || (length(value) > 0 && is.null(names(value))) ||
      any(names(value) == "") || any(duplicated(names(value)))) {
    stop("Value must be a named list with unique names")
  }
  if ( any(sapply(value, length) != vcount(graph)) ) {
    stop("Invalid attribute value length, must match number of vertices")
  }
  
  ## Trick to make R copy the graph
  newgraph <- graph
  attr(newgraph, "foo") <- NULL
  
  ## !!! Modifies the graph in place
  .Call("R_igraph_mybracket2_set", newgraph, 9L, 3L, value,
        PACKAGE="igraph")
  
  newgraph
}

get.edge.attribute <- function(graph, name, index=E(graph)) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  name <- as.character(name)
  index <- as.igraph.es(graph, index)
  myattr <- .Call("R_igraph_mybracket2", graph, 9L, 4L,
                  PACKAGE="igraph")[[name]]
  if (is.list(myattr) && length(index)==1) {
    myattr[[index]]
  } else {
    myattr[index]
  }
}

set.edge.attribute <- function(graph, name, index=E(graph), value) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  name <- as.character(name)
  index <- as.igraph.es(graph, index)
  ec <- ecount(graph)

  eattrs <- .Call("R_igraph_mybracket2", graph, 9L, 4L, PACKAGE="igraph")
  eattrs[[name]][index] <- value
  length(eattrs[[name]]) <- ec

  ## Trick to make R copy the graph
  newgraph <- graph
  attr(newgraph, "foo") <- NULL

  ## !!! Modifies the graph in place
  .Call("R_igraph_mybracket2_set", newgraph, 9L, 4L, eattrs, PACKAGE="igraph")

  newgraph
}

edge.attributes <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_mybracket2_copy", graph, 9L, 4L, PACKAGE="igraph")  
}

"edge.attributes<-" <- function(graph, value) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (!is.list(value) || (length(value) > 0 && is.null(names(value))) ||
      any(names(value) == "") || any(duplicated(names(value)))) {
    stop("Value must be a named list with unique names")
  }
  if ( any(sapply(value, length) != ecount(graph)) ) {
    stop("Invalid attribute value length, must match number of edges")
  }
  
  ## Trick to make R copy the graph
  newgraph <- graph
  attr(newgraph, "foo") <- NULL
  
  ## !!! Modifies the graph in place
  .Call("R_igraph_mybracket2_set", newgraph, 9L, 4L, value,
        PACKAGE="igraph")
  
  newgraph
}

list.graph.attributes <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  res <- names(.Call("R_igraph_mybracket2", graph, 9L, 2L, PACKAGE="igraph"))
  if (is.null(res)) { res <- character() }
  res
}

list.vertex.attributes <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  res <- names(.Call("R_igraph_mybracket2", graph, 9L, 3L, PACKAGE="igraph"))
                     
  if (is.null(res)) { res <- character() }
  res
}

list.edge.attributes <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  res <- names(.Call("R_igraph_mybracket2", graph, 9L, 4L, PACKAGE="igraph"))
  if (is.null(res)) { res <- character() }
  res
}

remove.graph.attribute <- function(graph, name) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  name <- as.character(name)
  if (!name %in% list.graph.attributes(graph)) {
    stop("No such graph attribute: ", name)
  }

  gattr <- .Call("R_igraph_mybracket2", graph, 9L, 2L, PACKAGE="igraph")
  gattr[[name]] <- NULL
  
  ## Trick to make R copy the graph
  newgraph <- graph
  attr(newgraph, "foo") <- NULL

  ## !!! Modifies the graph in place
  .Call("R_igraph_mybracket2_set", newgraph, 9L, 2L, gattr, PACKAGE="igraph")

  newgraph
}

remove.vertex.attribute <- function(graph, name) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  name <- as.character(name)
  if (!name %in% list.vertex.attributes(graph)) {
    stop("No such vertex attribute: ", name)
  }

  vattr <- .Call("R_igraph_mybracket2", graph, 9L, 3L, PACKAGE="igraph")
  vattr[[name]] <- NULL
  
  ## Trick to make R copy the graph
  newgraph <- graph
  attr(newgraph, "foo") <- NULL

  ## !!! Modifies the graph in place
  .Call("R_igraph_mybracket2_set", newgraph, 9L, 3L, vattr, PACKAGE="igraph")

  newgraph
}

remove.edge.attribute <- function(graph, name) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  name <- as.character(name)
  if (!name %in% list.edge.attributes(graph)) {
    stop("No such edge attribute: ", name)
  }

  eattr <- .Call("R_igraph_mybracket2", graph, 9L, 4L, PACKAGE="igraph")
  eattr[[name]] <- NULL
  
  ## Trick to make R copy the graph
  newgraph <- graph
  attr(newgraph, "foo") <- NULL

  ## !!! Modifies the graph in place
  .Call("R_igraph_mybracket2_set", newgraph, 9L, 4L, eattr, PACKAGE="igraph")

  newgraph
}

#############

is.named <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  "name" %in% list.vertex.attributes(graph)
}

is.weighted <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  "weight" %in% list.edge.attributes(graph)
}

is.bipartite <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  "type" %in% list.vertex.attributes(graph)
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
