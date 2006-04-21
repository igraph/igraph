
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

get.graph.attribute <- function(graph) {
  graph[[9]][[2]]
}

set.graph.attribute <- function(graph, name, value) {
  graph[[9]][[2]][[name]] <- value
  graph
}

get.vertex.attribute <- function(graph, name=NULL) {
  graph[[9]][[3]]
}

set.vertex.attribute <- function(graph, name, value) {
  graph[[9]][[3]][[name]] <- value
  graph
}

get.edge.attribute <- function(graph, name=NULL) {
  graph[[9]][[4]]
}

set.edge.attribute <- function(graph, name, value) {
  graph[[9]][[4]][[name]] <- value
  graph
}

g <- get.graph.attribute
v <- get.vertex.attribute
e <- get.edge.attribute

