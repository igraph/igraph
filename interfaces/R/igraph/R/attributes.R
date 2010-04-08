
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
  graph[[9]][[2]][[as.character(name)]]
}

set.graph.attribute <- function(graph, name, value) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  graph[[9]][[2]][[as.character(name)]] <- value
  graph
}

get.vertex.attribute <- function(graph, name, index=V(graph)) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  index <- as.igraph.vs(graph, index)
  name <- as.character(name)
  if (is.list(graph[[9]][[3]][[name]]) && length(index)==1) {
    graph[[9]][[3]][[as.character(name)]][[index+1]]
  } else {
    graph[[9]][[3]][[as.character(name)]][index+1]
  }
}

set.vertex.attribute <- function(graph, name, index=V(graph), value) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  index <- as.igraph.vs(graph, index)
  name <- as.character(name)
  vc <- vcount(graph)
##   if (length(index) %% length(value)) {
##     warning("number of items to replace is not a multiple of replacement length")
##   }
##   value <- rep(value, length.out=length(index))
  graph[[9]][[3]][[name]][index+1] <- value
  length(graph[[9]][[3]][[name]]) <- vc
  graph
}

get.edge.attribute <- function(graph, name, index=E(graph)) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  name <- as.character(name)
  if (is.list(graph[[9]][[4]][[name]]) && length(index)==1) {
    graph[[9]][[4]][[name]][[index+1]]
  } else {
    graph[[9]][[4]][[name]][index+1]
  }
}

set.edge.attribute <- function(graph, name, index=E(graph), value) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  name <- as.character(name)
  ec <- ecount(graph)
##   if (length(index) %% length(value)) {
##     warning("number of items to replace is not a multiple of replacement length")
##   }
##   value <- rep(value, length.out=length(index))
  graph[[9]][[4]][[name]][index+1] <- value
  length(graph[[9]][[4]][[name]]) <- ec
  graph
}

list.graph.attributes <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  res <- names(graph[[9]][[2]])
  if (is.null(res)) { res <- character() }
  res
}

list.vertex.attributes <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  res <- names(graph[[9]][[3]])
  if (is.null(res)) { res <- character() }
  res
}

list.edge.attributes <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  res <- names(graph[[9]][[4]])
  if (is.null(res)) { res <- character() }
  res
}

remove.graph.attribute <- function(graph, name) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  graph[[9]][[2]][[as.character(name)]] <- NULL
  graph
}

remove.vertex.attribute <- function(graph, name) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  graph[[9]][[3]][[as.character(name)]] <- NULL
  graph
}

remove.edge.attribute <- function(graph, name) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  graph[[9]][[4]][[as.character(name)]] <- NULL
  graph
}

#############

igraph.i.attribute.combination <- function(comb) {
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
