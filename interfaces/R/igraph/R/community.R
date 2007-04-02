
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

###################################################################
# Community structure
###################################################################

spinglass.community <- function(graph, weights=NULL, vertex=NULL, spins=25,
                                parupdate=FALSE, start.temp=1,
                                stop.temp=0.01, cool.fact=0.99,
                                update.rule="config", gamma=1.0) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(weights)) {
    if ("weight" %in% list.edge.attributes(graph)) {
      weights <- as.numeric(E(g)$weight)
    } else {
      weights <- as.numeric(rep(1, ecount(graph)))
    }
  }

  if (is.character(update.rule)) {
    update.rule <- switch(update.rule, "simple"=0, "random"=0, "config"=1)
  }

  if (is.null(vertex)) {    
    .Call("R_igraph_spinglass_community", graph, weights,
          as.numeric(spins), as.logical(parupdate), as.numeric(start.temp),
          as.numeric(stop.temp), as.numeric(cool.fact),
          as.numeric(update.rule), as.numeric(gamma),
          PACKAGE="igraph")
  } else {
    .Call("R_igraph_spinglass_my_community", graph, weights,
          as.numeric(vertex), as.numeric(spins), 
          as.numeric(update.rule), as.numeric(gamma),
          PACKAGE="igraph")
  }    
}

walktrap.community <- function(graph, weights=E(g)$weight, steps=4, merges=TRUE,
                               modularity=FALSE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  if (!is.null(weights)) {
    weight <- as.numeric(weight)
  }

  .Call("R_igraph_walktrap_community", graph, weights, as.numeric(steps),
        as.logical(merges), as.logical(modularity),
        PACKAGE="igraph")
}

community.eb <- function(graph, directed=TRUE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  all.nodes <- vcount(graph)
  all.edges <- ecount(graph)
  res <- matrix(0, nc=2, nr=all.edges);

  for (i in 1:all.edges) {

    print(i)
    
    max.eb <- which.max(edge.betweenness(graph, directed=directed))
    res[i,1] <- (max.eb-1) %%  all.nodes + 1
    res[i,2] <- (max.eb-1) %/% all.nodes + 1

    graph <- delete.edges(graph, res[i,])
  }
  
  res
}

community.cut <- function(graph, edges, after.removing) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (after.removing >= nrow(edges)) { 
    res <- graph.empty(directed=is.directed(graph))
  } else {
    res <- graph( t(edges[(after.removing+1):nrow(edges),]), 
		  directed=is.directed(graph) )
  }

  res
}

edge.type.matrix <- function(graph, types) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (vcount(graph) != length(types)) {
    stop("'graph' and/or 'types' invalid, they should be of the same length")
  }

  no.of.types <- max(types)
  res <- matrix(0, nr=no.of.types, nc=no.of.types)

  el <- get.edgelist(graph, names=FALSE)
  if (length(el) != 0) {
    for (i in 1:nrow(el)) {
      res[ types[el[i,1]], types[el[i,2]] ] <-
        res[ types[el[i,1]], types[el[i,2]] ] + 1
    }
    res <- res / sum(res)
  }
  
  res
}

modularity <- function(graph, types) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  etm <- edge.type.matrix(graph, types)

  res <- sum(diag(etm)) - sum(etm %*% etm)
  
  res
}

community.eb2 <- function(graph, directed=TRUE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  all.nodes <- vcount(graph)
  all.edges <- ecount(graph)
  cno <- length(clusters(graph)$csize)
  res <- numeric()

  for (i in 1:all.edges) {

    print(i)
    
    max.eb <- which.max(edge.betweenness(graph, directed=directed))
    from <- (max.eb-1) %%  all.nodes + 1
    to <- (max.eb-1) %/% all.nodes + 1

    res <- c(res, from, to)
    graph <- delete.edges(graph, c(from, to))
    cno2 <- length(clusters(graph)$csize)
    if (cno2 > cno) { break }
  }
  
  matrix(res, nc=2, byrow=TRUE)
}
