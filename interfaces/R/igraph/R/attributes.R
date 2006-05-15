
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
  graph[[9]][[2]][[as.character(name)]]
}

set.graph.attribute <- function(graph, name, value) {
  graph[[9]][[2]][[as.character(name)]] <- value
  graph
}

get.vertex.attribute <- function(graph, name, index=V(g)) {
  name <- as.character(name)
  if (is.list(graph[[9]][[3]][[name]]) && length(index)==1) {
    graph[[9]][[3]][[as.character(name)]][[index+1]]
  } else {
    graph[[9]][[3]][[as.character(name)]][index+1]
  }
}

set.vertex.attribute <- function(graph, name, index=V(g), value) {
  name <- as.character(name)
  graph[[9]][[3]][[name]][index+1] <- value
  if (length(graph[[9]][[3]][[name]]) != vcount(graph)) {
    graph[[9]][[3]][[name]] <-
      rep(graph[[9]][[3]][[name]], length.out=vcount(graph))
  }
  graph
}

get.edge.attribute <- function(graph, name, index=E(g)) {
  name <- as.character(name)
  if (is.list(graph[[9]][[4]][[name]]) && length(index)==1) {
    graph[[9]][[4]][[name]][[index+1]]
  } else {
    graph[[9]][[4]][[name]][index+1]
  }
}

set.edge.attribute <- function(graph, name, index=E(g), value) {
  name <- as.character(name)
  graph[[9]][[4]][[name]][index+1] <- value
  if (length(graph[[9]][[4]][[name]]) != ecount(graph)) {
    graph[[9]][[4]][[name]] <-
      rep(graph[[9]][[4]][[name]], length.out=ecount(graph))
  }
  graph
}

list.graph.attributes <- function(graph) {
  res <- names(graph[[9]][[2]])
  if (is.null(res)) { res <- character() }
  res
}

list.vertex.attributes <- function(graph) {
  res <- names(graph[[9]][[3]])
  if (is.null(res)) { res <- character() }
  res
}

list.edge.attributes <- function(graph) {
  res <- names(graph[[9]][[4]])
  if (is.null(res)) { res <- character() }
  res
}

remove.graph.attribute <- function(graph, name) {
  graph[[9]][[2]][[as.character(name)]] <- NULL
  graph
}

remove.vertex.attribute <- function(graph, name) {
  graph[[9]][[3]][[as.character(name)]] <- NULL
  graph
}

remove.edge.attribute <- function(graph, name) {
  graph[[9]][[4]][[as.character(name)]] <- NULL
  graph
}
