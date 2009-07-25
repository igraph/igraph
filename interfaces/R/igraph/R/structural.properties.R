
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

diameter <- function(graph, directed=TRUE, unconnected=TRUE, weights=NULL) {
  
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) {
    weights <- E(graph)$weight
  }
  if (!is.null(weights) && any(!is.na(weights))) {
    weights <- as.numeric(weights)
  } else {
    weights <- NULL
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_diameter", graph, as.logical(directed),
        as.logical(unconnected), weights,
        PACKAGE="igraph")
}

get.diameter <- function(graph, directed=TRUE, unconnected=TRUE,
                         weights=NULL) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) {
    weights <- E(graph)$weight
  }
  if (!is.null(weights) && any(!is.na(weights))) {
    weights <- as.numeric(weights)
  } else {
    weights <- NULL
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_get_diameter", graph, as.logical(directed),
        as.logical(unconnected), weights,
        PACKAGE="igraph")
}

farthest.nodes <- function(graph, directed=TRUE, unconnected=TRUE,
                           weights=NULL) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) {
    weights <- E(graph)$weight
  }
  if (!is.null(weights) && any(!is.na(weights))) {
    weights <- as.numeric(weights)
  } else {
    weights <- NULL
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_farthest_points", graph, as.logical(directed),
        as.logical(unconnected), weights,
        PACKAGE="igraph")
}       

average.path.length <- function(graph, directed=TRUE, unconnected=TRUE) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_average_path_length", graph, as.logical(directed),
        as.logical(unconnected),
        PACKAGE="igraph")
}

degree <- function(graph, v=V(graph),
                   mode=c("all", "out", "in", "total"), loops=TRUE){
  
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3, "total"=3)
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
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

shortest.paths <- function(graph, v=V(graph), to=V(graph),
                           mode=c("all", "out", "in"),
                           weights=NULL,
                           algorithm=c("automatic", "unweighted", "dijkstra",
                             "bellman-ford", "johnson")) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)  
  algorithm <- igraph.match.arg(algorithm)
  algorithm <- switch(algorithm, "automatic"=0, "unweighted"=1,
                      "dijkstra"=2, "bellman-ford"=3, "johnson"=4)
  
  if (is.null(weights)) {
    if ("weight" %in% list.edge.attributes(graph)) {
      weights <- as.numeric(E(graph)$weight)
    }
  } else {
    if (length(weights)==1 && is.na(weights)) {
      weights <- NULL
    } else {
      weights <- as.numeric(weights)
    }
  }

  if (! is.null(weights) && algorithm==1) {
    weights <- NULL
    warning("Unweighted algorithm chosen, weights ignored")
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_shortest_paths", graph, as.igraph.vs(v),
        as.igraph.vs(to), as.numeric(mode), weights, as.numeric(algorithm),
        PACKAGE="igraph")
}

get.shortest.paths <- function(graph, from, to=V(graph),
                               mode=c("all", "out", "in"),
                               weights=NULL) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)

  if (is.null(weights)) {
    if ("weight" %in% list.edge.attributes(graph)) {
      weights <- as.numeric(E(graph)$weight)
    }
  } else {
    if (length(weights)==1 && is.na(weights)) {
      weights <- NULL
    } else {
      weights <- as.numeric(weights)
    }
  }
  
  to <- as.igraph.vs(to)
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_get_shortest_paths", graph,
        as.numeric(from), to, as.numeric(mode), as.numeric(length(to)),
        weights, PACKAGE="igraph")
}

get.all.shortest.paths <- function(graph, from,
                                   to=V(graph),
                                   mode=c("all", "out", "in")) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_get_all_shortest_paths", graph,
        as.numeric(from), as.igraph.vs(to), as.numeric(mode),
        PACKAGE="igraph")
}

subcomponent <- function(graph, v, mode=c("all", "out", "in")) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_subcomponent", graph, as.igraph.vs(v), as.numeric(mode),
        PACKAGE="igraph")
}

subgraph <- function(graph, v) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_subgraph", graph, as.igraph.vs(v),
        PACKAGE="igraph")
}

simplify <- function(graph, remove.multiple=TRUE,
                     remove.loops=TRUE) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
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

betweenness <- function(graph, v=V(graph), directed=TRUE, weights=NULL,
                        verbose=igraph.par("verbose")) {
  
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) {
    weights <- E(graph)$weight
  }
  if (!is.null(weights) && any(!is.na(weights))) {
    weights <- as.numeric(weights)
  } else {
    weights <- NULL
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_betweenness", graph, as.igraph.vs(v),
        as.logical(directed), weights, as.logical(verbose),
        PACKAGE="igraph")
}

transitivity <- function(graph, type=c("undirected", "global", "globalundirected",
                                  "localundirected", "local", "average",
                                  "localaverage", "localaverageundirected",
                                  "barrat", "weighted"),
                         vids=NULL, weights=NULL, isolates=c("NaN", "zero")) {
  
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  type <- igraph.match.arg(type)
  type <- switch(type, "undirected"=0, "global"=0, "globalundirected"=0,
                 "localundirected"=1, "local"=1, "average"=2,
                 "localaverage"=2, "localaverageundirected"=2, "barrat"=3,
                 "weighted"=3)

  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) {
    weights <- E(graph)$weight
  }
  if (!is.null(weights) && any(!is.na(weights))) {
    weights <- as.numeric(weights)
  } else {
    weights <- NULL
  }

  isolates <- igraph.match.arg(isolates)
  isolates <- as.double(switch(isolates, "nan"=0, "zero"=1))

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  if (type==0) {
    .Call("R_igraph_transitivity_undirected", graph, isolates,
          PACKAGE="igraph")
  } else if (type==1) {
    if (is.null(vids)) {
      .Call("R_igraph_transitivity_local_undirected_all", graph, isolates,
            PACKAGE="igraph")
    } else {
      vids <- as.igraph.vs(vids)
      .Call("R_igraph_transitivity_local_undirected", graph, as.numeric(vids),
            isolates, PACKAGE="igraph")
    }
  } else if (type==2) {
    .Call("R_igraph_transitivity_avglocal_undirected", graph, isolates,
          PACKAGE="igraph")
  } else if (type==3) {
    vids <- as.igraph.vs(vids)
    if (is.null(weights)) {
      .Call("R_igraph_transitivity_local_undirected", graph, as.numeric(vids),
            isolates, PACKAGE="igraph")
    } else { 
      .Call("R_igraph_transitivity_barrat", graph, as.numeric(vids), weights,
            isolates, PACKAGE="igraph")
    }
  }
}

## Generated by stimulus now
## graph.laplacian <- function(graph, normalized=FALSE) {

##   if (!is.igraph(graph)) {
##     stop("Not a graph object")
##   }
  
##   on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
##   .Call("R_igraph_laplacian", graph, as.logical(normalized),
##         PACKAGE="igraph")
## }
  
## OLD implementation
## graph.laplacian <- function(graph, normalized=FALSE) {

##   if (!is.igraph(graph)) {
##     stop("Not a graph object")
##   }
##   if (is.directed(graph)) {
##     warning("Laplacian of a directed graph???")
##   }

##   M <- get.adjacency(graph)
##   if (!normalized) {
##     M <- structure(ifelse(M>0, -1, 0), dim=dim(M))
##     diag(M) <- degree(graph)
##   } else {
##     deg <- degree(graph)
##     deg <- outer(deg, deg, "*")
##     M <- structure(ifelse(M>0, -1/deg, 0))
##     diag(M) <- 1
##   }
  
##   M
## }

## Structural holes a'la Burt, code contributed by
## Jeroen Bruggeman

## constraint.orig <- function(graph, nodes=V(graph), attr=NULL) {

##   if (!is.igraph(graph)) {
##     stop("Not a graph object")
##   }

##   idx <- degree(graph) != 0
##   A <- get.adjacency(graph, attr=attr)
##   A <- A[idx, idx]
##   n <- sum(idx)
  
##   one <- c(rep(1,n))
##   CZ <- A + t(A)
##   cs <- CZ %*% one                      # degree of vertices
##   ics <- 1/cs
##   CS <- ics %*% t(one)                  # 1/degree of vertices
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

## Newest implementation, hopefully correct, there is a C implementation
## now so we don't need this

## constraint.old <- function(graph, nodes=V(graph)) {

##   if (!is.igraph(graph)) {
##     stop("Not a graph object")
##   }

##   nodes <- as.numeric(nodes)
##   res <- numeric(length(nodes))
##   deg <- degree(graph, mode="all", loops=FALSE)

##   not <- function(i, v) v[ v!=i ]

##   for (a in seq(along=nodes)) {
##     i <- nodes[a]
    
##     first <- not(i, neighbors(graph, i, mode="all"))
##     first <- unique(first)
##     for (b in seq(along=first)) {
##       j <- first[b]

##       ## cj is the contribution of j
##       cj <- are.connected(graph, i, j)      / deg[i+1]
##       cj <- cj + are.connected(graph, j, i) / deg[i+1]

##       second <- not(i, not(j, neighbors(graph, j, mode="all")))
##       for (c in seq(along=second)) {
##         q <- second[c]
##         cj <- cj + are.connected(graph, i, q) / deg[q+1] / deg[i+1]
##         cj <- cj + are.connected(graph, q, i) / deg[q+1] / deg[i+1]
##       }
                            
##       ## Ok, we have the total contribution of j
##       res[a] <- res[a] + cj*cj
##     }
##   }

##   if (!is.directed(graph)) {
##     res <- res/4
##   }
##   res
## }

constraint <- function(graph, nodes=V(graph), weights=NULL) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(weights)) {
    if ("weight" %in% list.edge.attributes(graph)) {
      weights <- E(graph)$weight
    }
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_constraint", graph, as.igraph.vs(nodes),
        as.numeric(weights),        
        PACKAGE="igraph")
}

reciprocity <- function(graph, ignore.loops=TRUE) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_reciprocity", graph, as.logical(ignore.loops),
        PACKAGE="igraph")
}

rewire <- function(graph, mode="simple", niter=100) {
  
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "simple"=0)
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_rewire", graph, as.numeric(niter), as.numeric(mode),
        PACKAGE="igraph")
}

bonpow.dense <- function(graph, nodes=V(graph),
                         loops=FALSE, exponent=1,
                         rescale=FALSE, tol=1e-7){

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }  
  
  d <- get.adjacency(graph)
  if (!loops) {
    diag(d) <- 0
  }
  n <- vcount(graph)
  id <- matrix(0,nrow=n,ncol=n)
  diag(id) <- 1

#  ev <- apply(solve(id-exponent*d,tol=tol)%*%d,1,sum)
  ev <- solve(id-exponent*d, tol=tol) %*% apply(d,1,sum)
  if(rescale) {
    ev <- ev/sum(ev)
  } else {
    ev <- ev*sqrt(n/sum((ev)^2))
  } 
  ev[as.numeric(nodes)+1]
}

bonpow.sparse <- function(graph, nodes=V(graph), loops=FALSE,
                          exponent=1, rescale=FALSE, tol=1e-07) {

  ## remove loops if requested
  if (!loops) {
    graph <- simplify(graph, remove.multiple=FALSE, remove.loops=TRUE)
  }

  vg <- vcount(graph)
  
  ## sparse adjacency matrix
  d <- get.adjacency(graph, sparse=TRUE)

  ## sparse identity matrix
  id <- Diagonal(vg)

  ## solve it
  ev <- solve(id - exponent * d, degree(graph, mode="out"), tol=tol)

  if (rescale) {
    ev <- ev/sum(ev)
  } else {
    ev <- ev * sqrt(vcount(graph)/sum((ev)^2))
  }

  ev[as.numeric(nodes) + 1]
}

bonpow <- function(graph, nodes=V(graph),
                   loops=FALSE, exponent=1,
                   rescale=FALSE, tol=1e-7, sparse=TRUE){
  if (sparse) {
    if (require(Matrix)) {
      return(bonpow.sparse(graph, nodes, loops, exponent, rescale, tol))
    }
  } 

  bonpow.dense(graph, nodes, loops, exponent, rescale, tol)
}

alpha.centrality.dense <- function(graph, nodes=V(graph), alpha=1,
                                   loops=FALSE, exo=1, weights=NULL,
                                   tol=1e-7) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  exo <- rep(exo, length=vcount(graph))
  exo <- matrix(exo, nc=1)

  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) {
    ## weights == NULL and there is a "weight" edge attribute
    attr <- "weight"
  } else if (is.null(weights)) {
    ## weights == NULL, but there is no "weight" edge attribute
    attr <- NULL
  } else if (any(!is.na(weights))) {
    ## weights != NULL and weights != rep(NA, x)
    graph <- set.edge.attribute(graph, "weight", value=as.numeric(weights))
    attr <- "weight"
  } else {
    ## weights != NULL, but weights == rep(NA, x)
    attr <- NULL
  }

  d <- t(get.adjacency(graph, attr=attr))
  if (!loops) {
    diag(d) <- 0
  }
  n <- vcount(graph)
  id <- matrix(0, nrow=n, ncol=n)
  diag(id) <- 1
  
  ev <- solve(id-alpha*d, tol=tol) %*% exo
  ev[as.numeric(nodes)+1]
}

alpha.centrality.sparse <- function(graph, nodes=V(graph), alpha=1,
                                   loops=FALSE, exo=1, weights=NULL,
                                   tol=1e-7) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  vc <- vcount(graph)

  if (!loops) {
    graph <- simplify(graph, remove.multiple=FALSE, remove.loops=TRUE)
  }  

  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) {
    ## weights == NULL and there is a "weight" edge attribute
    weights <- E(graph)$weight
  } else if (is.null(weights)) {
    ## weights == NULL, but there is no "weight" edge attribute
    weights <- rep(1, ecount(graph))
  } else if (any(!is.na(weights))) {
    weights <- as.numeric(weights)
    ## weights != NULL, but weights == rep(NA, x)
    weights <- rep(1, ecount(graph))
  }
  
  el <- get.edgelist(graph, names=FALSE)+1
  M <- spMatrix(vc, vc, i=el[,2], j=el[,1], x=weights)
  M <- as(M, "dgCMatrix")
  
  ## Create an identity matrix
  M2 <- spMatrix(vc, vc, i=1:vc, j=1:vc, x=rep(1, vc))
  M2 <- as(M2, "dgCMatrix")

  ## exo
  exo <- cbind(rep(exo, length=vc))

  ## Solve the equation
  M3 <- M2-alpha*M
  r <- solve(M3, tol=tol, exo)
  
  r[ as.numeric(nodes)+1]
}

alpha.centrality <- function(graph, nodes=V(graph), alpha=1,
                             loops=FALSE, exo=1, weights=NULL,
                             tol=1e-7, sparse=TRUE) {
  if (sparse) {
    if (require(Matrix)) {
      return(alpha.centrality.sparse(graph, nodes, alpha, loops,
                                     exo, weights, tol))
    }
  } 
  alpha.centrality.dense(graph, nodes, alpha, loops, exo, weights, tol)  
}


graph.density <- function(graph, loops=FALSE) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }  
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_density", graph, as.logical(loops),
        PACKAGE="igraph")
}

neighborhood.size <- function(graph, order, nodes=V(graph),
                              mode=c("all", "out", "in")) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_neighborhood_size", graph, 
        as.igraph.vs(nodes), as.numeric(order), as.numeric(mode),
        PACKAGE="igraph")
}

neighborhood <- function(graph, order, nodes=V(graph), mode=c("all", "out", "in")) {
  
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_neighborhood", graph, 
        as.igraph.vs(nodes), as.numeric(order), as.numeric(mode),
        PACKAGE="igraph")
}

graph.neighborhood <- function(graph, order, nodes=V(graph),
                               mode=c("all", "out", "in")) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_neighborhood_graphs", graph, 
        as.igraph.vs(nodes), as.numeric(order), as.numeric(mode),
        PACKAGE="igraph")
}

graph.coreness <- function(graph, mode=c("all", "out", "in")) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_coreness", graph, as.numeric(mode),
        PACKAGE="igraph")
}

topological.sort <- function(graph, mode=c("out", "all", "in")) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_topological_sorting", graph, as.numeric(mode),
        PACKAGE="igraph")
}

girth <- function(graph, circle=TRUE) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_girth", graph, as.logical(circle),
        PACKAGE="igraph")
}

is.loop <- function(graph, eids=E(graph)) {

  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_is_loop", graph, as.igraph.es(eids),
        PACKAGE="igraph")
}

is.multiple <- function(graph, eids=E(graph)) {

  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_is_multiple", graph, as.igraph.es(eids),
        PACKAGE="igraph")
}

count.multiple <- function(graph, eids=E(graph)) {

  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_count_multiple", graph, as.igraph.es(eids),
        PACKAGE="igraph")
}

graph.bfs <- function(graph, root, neimode=c("out", "in", "all", "total"),
                      order=TRUE, rank=FALSE, father=FALSE,
                      pred=FALSE, succ=FALSE, dist=FALSE,
                      callback=NULL, extra=NULL, rho=parent.frame()) {

  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }

  root <- as.numeric(root)
  neimode <- switch(igraph.match.arg(neimode),
                    "out"=1, "in"=2, "all"=3, "total"=3)
  if (!is.null(callback)) { callback <- as.function(callback) }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_bfs", graph, root, neimode,
        as.logical(order), as.logical(rank), as.logical(father),
        as.logical(pred), as.logical(succ), as.logical(dist),
        callback, extra, rho,
        PACKAGE="igraph")
  
}

graph.dfs <- function(graph, root, neimode=c("out", "in", "all", "total"),
                      order=TRUE, order.out=FALSE, father=FALSE, dist=FALSE,
                      in.callback=NULL, out.callback=NULL, extra=NULL,
                      rho=parent.frame()) {

  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }

  root <- as.numeric(root)
  neimode <- switch(igraph.match.arg(neimode),
                    "out"=1, "in"=2, "all"=3, "total"=3)
  if (!is.null(in.callback)) { in.callback <- as.function(in.callback) }
  if (!is.null(out.callback)) { out.callback <- as.function(out.callback) }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_dfs", graph, root, neimode,
        as.logical(order), as.logical(order.out), as.logical(father),
        as.logical(dist), in.callback, out.callback, extra, rho,
        PACKAGE="igraph")
  
}

edge.betweenness <- function(graph, e=E(graph),
                             directed=TRUE, weights=NULL) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  e <- as.igraph.es(e)
  directed <- as.logical(directed)
  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) { 
  weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
  weights <- as.numeric(weights) 
  } else { 
  weights <- NULL 
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_edge_betweenness", graph, directed, weights,
        PACKAGE="igraph")
  res[as.numeric(e)+1]
}

edge.betweenness.estimate <- function(graph, e=E(graph),
                                      directed=TRUE, cutoff, weights=NULL) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  e <- as.igraph.es(e)
  directed <- as.logical(directed)
  cutoff <- as.numeric(cutoff)
  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) { 
  weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
  weights <- as.numeric(weights) 
  } else { 
  weights <- NULL 
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_edge_betweenness_estimate", graph, directed, cutoff, weights,
        PACKAGE="igraph")
  res[as.numeric(e)+1]
}
