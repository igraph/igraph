
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
# Structural properties
###################################################################

diameter <- function(graph, directed=TRUE, unconnected=TRUE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_diameter", graph, as.logical(directed),
        as.logical(unconnected),
        PACKAGE="igraph")
}

get.diameter <- function(graph, directed=TRUE, unconnected=TRUE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_get_diameter", graph, as.logical(directed),
        as.logical(unconnected),
        PACKAGE="igraph")
}

farthest.nodes <- function(graph, directed=TRUE, unconnected=TRUE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_farthest_points", graph, as.logical(directed),
        as.logical(unconnected),
        PACKAGE="igraph")
}       

average.path.length <- function(graph, directed=TRUE, unconnected=TRUE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_average_path_length", graph, as.logical(directed),
        as.logical(unconnected),
        PACKAGE="igraph")
}

degree <- function(graph, v=V(graph),
                   mode="total", loops=TRUE){
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (is.character(mode)) {
    mode <- switch(mode, "out"=1, "in"=2, "all"=3, "total"=3)
  }
  
  .Call("R_igraph_degree", graph, as.igraph.vs(v), as.numeric(mode),
        as.logical(loops), PACKAGE="igraph")
}
  
degree.distribution <- function(graph, cumulative=FALSE, ...) {
  
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  cs <- degree(graph, ...)
  hi <- hist(cs, -1:max(cs), plot=FALSE)$intensities
  if (!cumulative) {
    res <- hi
  } else {
    res <- rev(cumsum(rev(hi)))
  }
  
  res
}

closeness <- function(graph, v=V(graph), mode="all") {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (is.character(mode)) {
    mode <- switch(mode, "out"=1, "in"=2, "all"=3)
  }
  
  .Call("R_igraph_closeness", graph, as.igraph.vs(v), as.numeric(mode),
        PACKAGE="igraph")
}

shortest.paths <- function(graph, v=V(graph), mode="all") {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (is.character(mode)) {
    mode <- switch(mode, "out"=1, "in"=2, "all"=3)
  }

  .Call("R_igraph_shortest_paths", graph, as.igraph.vs(v),
        as.numeric(mode),
        PACKAGE="igraph")
}

get.shortest.paths <- function(graph, from, to=V(graph),
                               mode="all") {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (is.character(mode)) {
    mode <- switch(mode, "out"=1, "in"=2, "all"=3)
  }

  to <- as.igraph.vs(to)
  .Call("R_igraph_get_shortest_paths", graph,
        as.numeric(from), to, as.numeric(mode), as.numeric(length(to)),
        PACKAGE="igraph")
}

get.all.shortest.paths <- function(graph, from,
                                   to=V(graph),
                                   mode="all") {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (is.character(mode)) {
    mode <- switch(mode, "out"=1, "in"=2, "all"=3)
  }

  .Call("R_igraph_get_all_shortest_paths", graph,
        as.numeric(from), as.igraph.vs(to), as.numeric(mode),
        PACKAGE="igraph")
}

subcomponent <- function(graph, v, mode="all") {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (is.character(mode)) {
    mode <- switch(mode, "out"=1, "in"=2, "all"=3)
  }

  .Call("R_igraph_subcomponent", graph, as.igraph.vs(v), as.numeric(mode),
        PACKAGE="igraph")
}

subgraph <- function(graph, v) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_subgraph", graph, as.igraph.vs(v),
        PACKAGE="igraph")
}

simplify <- function(graph, remove.multiple=TRUE,
                     remove.loops=TRUE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_simplify", graph, as.logical(remove.multiple),
        as.logical(remove.loops), PACKAGE="igraph")

##   res <- graph
##   vc <- vcount(res)
##   if (remove.loops && vc > 0) {
##     remove <- numeric()
##     for (i in 0:(vc-1)) {
##       neis <- neighbors(graph, i, "out")
##       loops <- sum(neis==i)
##       if (is.directed(graph)) { loops <- loops*2 }
##       remove <- c(remove, rep(i, loops))
##     }
##     res <- delete.edges(res, remove)  
##   }
##   if (remove.multiple) {
##     remove <- numeric()
##     for (i in 0:(vc-1)) {
##       neis <- neighbors(graph, i, "out")
##       dup <- neis[ duplicated(neis) & neis > i ]
##       l <- sum(neis==i)
##       if (l>2) { dup <- c(dup, rep(i, l/4)) }
##       remove <- c(remove, as.numeric(t(matrix(c(rep(i,length(dup)),
##                                                 dup), nc=2))))
##     }
##     res <- delete.edges(res, remove)  
##   }  
  
##   res
}

betweenness <- function(graph, v=V(graph), directed=TRUE) {
  
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_betweenness", graph, as.igraph.vs(v),
        as.logical(directed),
        PACKAGE="igraph")
}

edge.betweenness <- function(graph, e=E(graph), directed=TRUE) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  .Call("R_igraph_edge_betweenness", graph, as.logical(directed),
        PACKAGE="igraph")[ as.numeric(e)+1 ]  
}

transitivity <- function(graph, type="undirected") {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (is.character(type)) {
    type <- switch(type, "undirected"=0)
  }

  .Call("R_igraph_transitivity", graph, as.numeric(type),
        PACKAGE="igraph")
}

graph.laplacian <- function(graph, normalized=FALSE) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (is.directed(graph)) {
    warning("Laplacian of a directed graph???")
  }

  M <- get.adjacency(graph)
  if (!normalized) {
    M <- structure(ifelse(M>0, -1, 0), dim=dim(M))
    diag(M) <- degree(graph)
  } else {
    deg <- degree(graph)
    deg <- outer(deg, deg, "*")
    M <- structure(ifelse(M>0, -1/deg, 0))
    diag(M) <- 1
  }
  
  M
}

## Structural holes a'la Burt, code contributed by
## Jeroen Bruggeman

## constraint <- function(graph, nodes=V(graph)) {

##   if (!is.igraph(graph)) {
##     stop("Not a graph object")
##   }

##   idx <- degree(graph) != 0
##   A <- get.adjacency(graph)
##   A <- A[idx, idx]
##   n <- sum(idx)
  
##   one <- c(rep(1,n))
##   CZ <- A + t(A)
##   cs <- CZ %*% one
##   ics <- 1/cs
##   CS <- ics %*% t(one)
##   P <- CZ * CS  #intermediate result: proportionate tie strengths
##   PSQ <- P%*%P #sum paths of length two
##   P.bi <- as.numeric(P>0)  #exclude paths to non-contacts (& reflexive):
##   PC <- (P + (PSQ*P.bi))^2  #dyadic constraint
##   ci <- PC %*% one      #overall constraint
##   dim(ci) <- NULL

##   ci2 <- numeric(vcount(graph))
##   ci2[idx] <- ci
##   ci2[!idx] <- NaN
##   ci2[nodes+1]
## }

## Second implementation

## constraint <- function(graph, nodes=V(graph)) {

##   if (!is.igraph(graph)) {
##     stop("Not a graph object")
##   }

##   nodes <-as.numeric(nodes)
  
##   weights <- 1/degree(graph)
##   res <- numeric(length(nodes))
##   for (i in seq(along=nodes)) {
##     res[i] <- res[i] + weights[nodes[i]+1]
##     first <- neighbors(graph, nodes[i], mode="all")
##     first <- first [ first != nodes[i] ]
##     for (j in seq(along=first)) {
##       second <- neighbors(graph, first[j], mode="all")
##       second <- second [ second %in% first ]
##       res[i] <- res[i] + sum(weights[first[j]+1] * weights[second+1])
##     }
##   }
##   res [ res == Inf ] <- NaN
##   res
## }

constraint <- function(graph, nodes=V(graph)) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  .Call("R_igraph_constraint", graph, as.igraph.vs(nodes),
        PACKAGE="igraph")
}
