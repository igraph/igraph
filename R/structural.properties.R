
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
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
###################################################################

###################################################################
# Structural properties
###################################################################

diameter <- function(graph, directed=TRUE, unconnected=TRUE) {

  res <- .Call("REST_diameter", igraph.c.interface,
               graph, as.logical(directed), as.logical(unconnected),
               PACKAGE="igraph")
  res
}

degree <- function(graph, v=1:vcount(graph), mode="total", loops=TRUE) {

  res <- numeric(length(v))
  if (!is.directed(graph) || mode %in% c("out", "total")) {
    for (i in seq(along=v)) {
      tmp <- neighbors(graph, v[i], "out")
      if (!loops)
        { res[i] <- sum(tmp != v[i]) }
      else
        { res[i] <- length(tmp) }
    }
  }
  if (is.directed(graph) && mode %in% c("in", "total")) {
    for (i in seq(along=v)) {
      tmp <- neighbors(graph, v[i], "in")
      if (!loops)
        { res[i] <- res[i] + sum(tmp != v[i]) }
      else
        { res[i] <- res[i]+ length(tmp) }
    }
  }

  res
}

degree.distribution <- function(graph, cumulative=FALSE, ...) {
  
  cs <- degree(graph, ...)
  hi <- hist(cs, -1:max(cs), plot=FALSE)$intensities
  if (!cumulative) {
    res <- hi
  } else {
    res <- rev(cumsum(rev(hi)))
  }
  
  res
}

closeness <- function(graph, v=1:vcount(graph), mode="all") {

  res <- .Call("REST_closeness", igraph.c.interface,
               graph, as.numeric(v), as.character(mode),
               PACKAGE="igraph")
  res
}

shortest.paths <- function(graph, v=1:vcount(graph), mode="all") {

  res <- .Call("REST_shortest_paths", igraph.c.interface, graph,
               as.double(v), mode,
               PACKAGE="igraph")
  
  res
}

subcomponent <- function(graph, v, mode="all") {

  res <- .Call("REST_subcomponent", igraph.c.interface, graph, v, mode,
               PACKAGE="igraph")
  
  res
}

subgraph <- function(graph, v) {
  res <- delete.vertices(graph, (1:vcount(graph))[-v])
  res
}

simplify <- function(graph, remove.loops=TRUE,
                     remove.multiple=TRUE) {

  res <- graph
  vc <- vcount(res)
  if (remove.loops && vc > 0) {
    remove <- numeric()
    for (i in 1:vc) {
      neis <- neighbors(graph, i, "out")
      loops <- sum(neis==i)
      if (is.directed(graph)) { loops <- loops*2 }
      remove <- c(remove, rep(i, loops))
    }
    res <- delete.edges(res, remove)  
  }
  if (remove.multiple) {
    remove <- numeric()
    for (i in 1:vc) {
      neis <- neighbors(graph, i, "out")
      dup <- neis[ duplicated(neis) ]
      remove <- c(remove, as.numeric(t(matrix(c(rep(i,length(dup)),
                                                dup), nc=2))))
    }
    res <- delete.edges(res, remove)  
  }  
  
  res
}

betweenness <- function(graph, v=1:vcount(graph), directed=TRUE) {

  res <- .Call("REST_betweenness", igraph.c.interface,
               graph, directed, PACKAGE="igraph")

  res <- res[v]
  
  if (!directed || !is.directed(graph)) {
    res <- res/2
  }

  res  
}
