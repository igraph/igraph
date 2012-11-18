
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
  res <- .Call("R_igraph_get_diameter", graph, as.logical(directed),
               as.logical(unconnected), weights,
               PACKAGE="igraph")
  res + 1
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
  res <- .Call("R_igraph_farthest_points", graph, as.logical(directed),
               as.logical(unconnected), weights,
               PACKAGE="igraph")
  res[1:2] <- res[1:2] + 1
  res
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
                   mode=c("all", "out", "in", "total"), loops=TRUE,
                   normalized=FALSE){
  
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  v <- as.igraph.vs(graph, v)
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3, "total"=3)
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_degree", graph, v-1,
               as.numeric(mode), as.logical(loops), PACKAGE="igraph")
  if (normalized) { res <- res / (vcount(graph)-1) }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[v]
  }
  res
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
  v <- as.igraph.vs(graph, v)
  to <- as.igraph.vs(graph, to)
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
  res <- .Call("R_igraph_shortest_paths", graph, v-1, to-1,
               as.numeric(mode), weights, as.numeric(algorithm),
               PACKAGE="igraph")
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    rownames(res) <- V(graph)$name[v]
    colnames(res) <- V(graph)$name[to]
  }
  res
}

get.shortest.paths <- function(graph, from, to=V(graph),
                               mode=c("out", "all", "in"),
                               weights=NULL,
                               output=c("vpath", "epath", "both"),
                               predecessors=FALSE, inbound.edges=FALSE) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)
  output <- igraph.match.arg(output)
  output <- switch(output, "vpath"=0, "epath"=1, "both"=2)

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
  
  to <- as.igraph.vs(graph, to)-1
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_get_shortest_paths", graph,
               as.igraph.vs(graph, from)-1, to, as.numeric(mode),
               as.numeric(length(to)), weights, as.numeric(output),
               as.logical(predecessors), as.logical(inbound.edges), 
               PACKAGE="igraph")

  if (!is.null(res$vpath)) {
    res$vpath <- lapply(res$vpath, function(x) x+1)
  }
  if (!is.null(res$epath)) {
    res$epath <- lapply(res$epath, function(x) x+1)
  }
  if (!is.null(res$predecessors)) {
    res$predecessors <- res$predecessors + 1
  }
  if (!is.null(res$inbound_edges)) {
    res$inbound_edges <- res$inbound_edges + 1
  }

  res
}

get.all.shortest.paths <- function(graph, from,
                                   to=V(graph),
                                   mode=c("out", "all", "in"),
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

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  if (is.null(weights)) {
    res <- .Call("R_igraph_get_all_shortest_paths", graph,
                 as.igraph.vs(graph, from)-1, as.igraph.vs(graph, to)-1,
                 as.numeric(mode), PACKAGE="igraph")
  } else {
    res <- .Call("R_igraph_get_all_shortest_paths_dijkstra", graph, 
                 as.igraph.vs(graph, from)-1, as.igraph.vs(graph, to)-1,
                 weights, as.numeric(mode), PACKAGE="igraph")
  }       
  res
}

subcomponent <- function(graph, v, mode=c("all", "out", "in")) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_subcomponent", graph, as.igraph.vs(graph, v)-1,
               as.numeric(mode),
               PACKAGE="igraph")
  res+1
}

subgraph <- function(graph, v) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_subgraph", graph, as.igraph.vs(graph, v)-1,
        PACKAGE="igraph")
}

betweenness <- function(graph, v=V(graph), directed=TRUE, weights=NULL,
                        nobigint=TRUE, normalized=FALSE) {
  
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  v <- as.igraph.vs(graph, v)
  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) {
    weights <- E(graph)$weight
  }
  if (!is.null(weights) && any(!is.na(weights))) {
    weights <- as.numeric(weights)
  } else {
    weights <- NULL
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_betweenness", graph, v-1,
               as.logical(directed), weights, as.logical(nobigint),
               PACKAGE="igraph")
  if (normalized) {
    vc <- vcount(graph)
    res <- 2*res / ( vc*vc-3*vc+2)
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[v]
  }
  res
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
      vids <- as.igraph.vs(graph, vids)-1
      .Call("R_igraph_transitivity_local_undirected", graph, vids,
            isolates, PACKAGE="igraph")
    }
  } else if (type==2) {
    .Call("R_igraph_transitivity_avglocal_undirected", graph, isolates,
          PACKAGE="igraph")
  } else if (type==3) {
    if (is.null(vids)) { vids <- V(graph) }
    vids <- as.igraph.vs(graph, vids)-1
    if (is.null(weights)) {
      .Call("R_igraph_transitivity_local_undirected", graph, vids,
            isolates, PACKAGE="igraph")
    } else { 
      .Call("R_igraph_transitivity_barrat", graph, vids, weights,
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
  nodes <- as.igraph.vs(graph, nodes)
  
  if (is.null(weights)) {
    if ("weight" %in% list.edge.attributes(graph)) {
      weights <- E(graph)$weight
    }
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_constraint", graph, nodes-1, as.numeric(weights),
               PACKAGE="igraph")
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[nodes]
  }
  res
}

reciprocity <- function(graph, ignore.loops=TRUE,
                        mode=c("default", "ratio")) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- switch(igraph.match.arg(mode), 'default'=0, 'ratio'=1)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_reciprocity", graph, as.logical(ignore.loops),
        as.numeric(mode), PACKAGE="igraph")
}

rewire <- function(graph, mode=c("simple", "loops"), niter=100) {
  
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "simple"=0, "loops"=1)
  
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
  ev[as.numeric(nodes)]
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

  ev[as.numeric(nodes)]
}

bonpow <- function(graph, nodes=V(graph),
                   loops=FALSE, exponent=1,
                   rescale=FALSE, tol=1e-7, sparse=TRUE){

  nodes <- as.igraph.vs(graph, nodes)
  if (sparse && require(Matrix)) {
    res <- bonpow.sparse(graph, nodes, loops, exponent, rescale, tol)
  }  else {
    res <- bonpow.dense(graph, nodes, loops, exponent, rescale, tol)
  }

  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- get.vertex.attribute(graph, "name", nodes)
  }
  
  res
}

alpha.centrality.dense <- function(graph, nodes=V(graph), alpha=1,
                                   loops=FALSE, exo=1, weights=NULL,
                                   tol=1e-7) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  exo <- rep(exo, length=vcount(graph))
  exo <- matrix(exo, ncol=1)

  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) {
    ## weights == NULL and there is a "weight" edge attribute
    attr <- "weight"
  } else if (is.null(weights)) {
    ## weights == NULL, but there is no "weight" edge attribute
    attr <- NULL
  } else if (is.character(weights) && length(weights)==1) {
    ## name of an edge attribute, nothing to do
    attr <- "weight"
  } else if (any(!is.na(weights))) {
    ## weights != NULL and weights != rep(NA, x)
    graph <- set.edge.attribute(graph, "weight", value=as.numeric(weights))
    attr <- "weight"
  } else {
    ## weights != NULL, but weights == rep(NA, x)
    attr <- NULL
  }

  d <- t(get.adjacency(graph, attr=attr, sparse=FALSE))
  if (!loops) {
    diag(d) <- 0
  }
  n <- vcount(graph)
  id <- matrix(0, nrow=n, ncol=n)
  diag(id) <- 1
  
  ev <- solve(id-alpha*d, tol=tol) %*% exo
  ev[as.numeric(nodes)]
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
  } else if (is.character(weights) && length(weights)==1) {
    weights <- get.edge.attribute(graph, weights)
  } else if (any(!is.na(weights))) {
    weights <- as.numeric(weights)
  } else {
    ## weights != NULL, but weights == rep(NA, x)
    weights <- rep(1, ecount(graph))
  } 
  
  el <- get.edgelist(graph, names=FALSE)
  M <- sparseMatrix(dims=c(vc, vc), i=el[,2], j=el[,1], x=weights)
  M <- as(M, "dgCMatrix")
  
  ## Create an identity matrix
  M2 <- sparseMatrix(dims=c(vc, vc), i=1:vc, j=1:vc, x=rep(1, vc))
  M2 <- as(M2, "dgCMatrix")

  ## exo
  exo <- cbind(rep(exo, length=vc))

  ## Solve the equation
  M3 <- M2-alpha*M
  r <- solve(M3, tol=tol, exo)
  
  r[ as.numeric(nodes)]
}

alpha.centrality <- function(graph, nodes=V(graph), alpha=1,
                             loops=FALSE, exo=1, weights=NULL,
                             tol=1e-7, sparse=TRUE) {

  nodes <- as.igraph.vs(graph, nodes)
  if (sparse && require(Matrix)) {
    res <- alpha.centrality.sparse(graph, nodes, alpha, loops,
                                   exo, weights, tol)
  } else {
    res <- alpha.centrality.dense(graph, nodes, alpha, loops,
                                  exo, weights, tol)
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- get.vertex.attribute(graph, "name", nodes)
  }
  res
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
        as.igraph.vs(graph, nodes)-1, as.numeric(order), as.numeric(mode),
        PACKAGE="igraph")
}

neighborhood <- function(graph, order, nodes=V(graph), mode=c("all", "out", "in")) {
  
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_neighborhood", graph, 
               as.igraph.vs(graph, nodes)-1, as.numeric(order),
               as.numeric(mode),
               PACKAGE="igraph")
  res <- lapply(res, function(x) x+1)
  res
}

graph.neighborhood <- function(graph, order, nodes=V(graph),
                               mode=c("all", "out", "in")) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_neighborhood_graphs", graph, 
               as.igraph.vs(graph, nodes)-1, as.numeric(order),
               as.numeric(mode),
               PACKAGE="igraph")
  res
}

graph.coreness <- function(graph, mode=c("all", "out", "in")) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_coreness", graph, as.numeric(mode),
               PACKAGE="igraph")
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- get.vertex.attribute(graph, "name")
  }
  res
}

topological.sort <- function(graph, mode=c("out", "all", "in")) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_topological_sorting", graph, as.numeric(mode),
               PACKAGE="igraph")
  res+1
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
  .Call("R_igraph_is_loop", graph, as.igraph.es(graph, eids)-1,
        PACKAGE="igraph")
}

is.multiple <- function(graph, eids=E(graph)) {

  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_is_multiple", graph, as.igraph.es(graph, eids)-1,
        PACKAGE="igraph")
}

count.multiple <- function(graph, eids=E(graph)) {

  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_count_multiple", graph, as.igraph.es(graph, eids)-1,
        PACKAGE="igraph")
}

graph.bfs <- function(graph, root, neimode=c("out", "in", "all", "total"),
                      unreachable=TRUE, restricted=NULL,
                      order=TRUE, rank=FALSE, father=FALSE,
                      pred=FALSE, succ=FALSE, dist=FALSE,
                      callback=NULL, extra=NULL, rho=parent.frame()) {

  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }

  if (length(root)==1) {
    root <- as.igraph.vs(graph, root)-1
    roots <- NULL    
  } else {
    root <- as.igraph.vs(graph, 0)      # ignored anyway
    roots <- as.igraph.vs(graph, root)
  }
  neimode <- switch(igraph.match.arg(neimode),
                    "out"=1, "in"=2, "all"=3, "total"=3)
  unreachable <- as.logical(unreachable)
  if (!is.null(restricted)) { restricted <- as.igraph.vs(graph, restricted) }
  if (!is.null(callback)) { callback <- as.function(callback) }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_bfs", graph, root, roots, neimode, unreachable,
               restricted,
               as.logical(order), as.logical(rank), as.logical(father),
               as.logical(pred), as.logical(succ), as.logical(dist),
               callback, extra, rho,
               PACKAGE="igraph")
  
  if (order)  res$order  <- res$order+1
  if (rank)   res$rank   <- res$rank+1
  if (father) res$father <- res$father+1
  if (pred)   res$pred   <- res$pred+1
  if (succ)   res$succ   <- res$succ+1
  res
}

graph.dfs <- function(graph, root, neimode=c("out", "in", "all", "total"),
                      unreachable=TRUE,
                      order=TRUE, order.out=FALSE, father=FALSE, dist=FALSE,
                      in.callback=NULL, out.callback=NULL, extra=NULL,
                      rho=parent.frame()) {

  if (!is.igraph(graph)) {
    stop("Not a graph object");
  }

  root <- as.igraph.vs(graph, root)-1
  neimode <- switch(igraph.match.arg(neimode),
                    "out"=1, "in"=2, "all"=3, "total"=3)
  unreachable <- as.logical(unreachable)
  if (!is.null(in.callback)) { in.callback <- as.function(in.callback) }
  if (!is.null(out.callback)) { out.callback <- as.function(out.callback) }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_dfs", graph, root, neimode, unreachable,
               as.logical(order), as.logical(order.out), as.logical(father),
               as.logical(dist), in.callback, out.callback, extra, rho,
               PACKAGE="igraph")
  
  if (order)     res$order     <- res$order+1
  if (order.out) res$order.out <- res$order.out+1
  if (father)    res$father    <- res$father+1
  res
}

edge.betweenness <- function(graph, e=E(graph),
                             directed=TRUE, weights=NULL) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  e <- as.igraph.es(graph, e)
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
  res[as.numeric(e)]
}

edge.betweenness.estimate <- function(graph, e=E(graph),
                                      directed=TRUE, cutoff, weights=NULL) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  e <- as.igraph.es(graph, e)
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
  res[as.numeric(e)]
}

clusters <- function(graph, mode=c("weak", "strong")) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  mode <- switch(igraph.match.arg(mode), "weak"=1, "strong"=2)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_clusters", graph, mode,
        PACKAGE="igraph")
  res$membership <- res$membership + 1
  res
}

unfold.tree <- function(graph, mode=c("all", "out", "in", "total"), roots) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  mode <- switch(igraph.match.arg(mode), "out"=1, "in"=2, "all"=3, "total"=3)
  roots <- as.igraph.vs(graph, roots)-1

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_unfold_tree", graph, mode, roots,
        PACKAGE="igraph")
  res
}

closeness <- function(graph, vids=V(graph),
                      mode=c("out", "in", "all", "total"), weights=NULL,
                      normalized=FALSE) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  vids <- as.igraph.vs(graph, vids)
  mode <- switch(igraph.match.arg(mode), "out"=1, "in"=2, "all"=3, "total"=3)
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
  res <- .Call("R_igraph_closeness", graph, vids-1, mode, weights,
               PACKAGE="igraph")
  if (!normalized) { res <- res / (vcount(graph)-1) }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }
  res
}

graph.laplacian <- function(graph, normalized=FALSE, weights=NULL,
                            sparse=getIgraphOpt("sparsematrices")) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  normalized <- as.logical(normalized)
  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) { 
    weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
    weights <- as.numeric(weights) 
  } else { 
    weights <- NULL 
  }
  sparse <- as.logical(sparse)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_laplacian", graph, normalized, weights, sparse,
               PACKAGE="igraph")
  if (sparse) {
    res <- igraph.i.spMatrix(res)
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    rownames(res) <- colnames(res) <- V(graph)$name
  }
  res
}

is.matching <- function(graph, matching, types=NULL) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  if (is.null(types) && "type" %in% list.vertex.attributes(graph)) { 
    types <- V(graph)$type 
  } 
  if (!is.null(types)) { 
    types <- as.logical(types) 
  }
  matching <- as.igraph.vs(graph, matching, na.ok=TRUE)-1
  matching[ is.na(matching) ] <- -1

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_is_matching", graph, types, matching,
        PACKAGE="igraph")

  res
}

is.maximal.matching <- function(graph, matching, types=NULL) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  if (is.null(types) && "type" %in% list.vertex.attributes(graph)) { 
    types <- V(graph)$type 
  } 
  if (!is.null(types)) { 
    types <- as.logical(types) 
  }
  matching <- as.igraph.vs(graph, matching, na.ok=TRUE)-1
  matching[ is.na(matching) ] <- -1

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_is_maximal_matching", graph, types, matching,
        PACKAGE="igraph")

  res
}

maximum.bipartite.matching <- function(graph, types=NULL, weights=NULL,
                                       eps=.Machine$double.eps) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  if (is.null(types) && "type" %in% list.vertex.attributes(graph)) { 
    types <- V(graph)$type 
  } 
  if (!is.null(types)) { 
    types <- as.logical(types) 
  }
  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) { 
    weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
    weights <- as.numeric(weights) 
  } else { 
    weights <- NULL 
  }
  eps <- as.numeric(eps)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_maximum_bipartite_matching", graph, types, weights,
               eps,
               PACKAGE="igraph")

  res$matching[ res$matching==0 ] <- NA
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    res$matching <- V(graph)$name[res$matching]
    names(res$matching) <- V(graph)$name
  }
  res
}
