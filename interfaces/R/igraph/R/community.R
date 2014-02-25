
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
# Community structure
###################################################################

membership <- function(communities) {
  if (!is.null(communities$membership)) {
    res <- communities$membership
  } else if (!is.null(communities$merges) &&
             !is.null(communities$modularity)) {
    res <- community.to.membership2(communities$merges, communities$vcount,
                                    which.max(communities$modularity))
  } else {
    stop("Cannot calculate community membership")
  }
  if (!is.null(communities$names)) {
    names(res) <- communities$names
  }
  res
}

print.communities <- function(x, ...) {
  cat("Graph community structure calculated with the",
      algorithm(x), "algorithm\n")
  if (algorithm(x)=="spinglass") {
    cat("Number of communities:", max(membership(x)), "\n")
    cat("Modularity:", modularity(x), "\n")
    cat("Membership vector:\n")
    print(membership(x))
  } else if (algorithm(x) %in% c("walktrap", "edge betweenness",
                                "fast greedy")) {
    cat("Number of communities (best split):", max(membership(x)), "\n")
    cat("Modularity (best split):", modularity(x), "\n")
    cat("Membership vector:\n")
    print(membership(x))
  } else if (algorithm(x) %in% c("leading eigenvector")) {
    cat("Number of communities (best split):", max(membership(x)), "\n")
    cat("Modularity (best split):", modularity(x), "\n")
    cat("Membership vector:\n")
    print(membership(x))
  } else if (algorithm(x) == "label propagation") {
    cat("Number of communities:", max(membership(x)), "\n")
    cat("Modularity:", modularity(x), "\n")
    cat("Membership vector:\n")
    print(membership(x))
  } else if (algorithm(x) == "multi level") {
    cat("Number of communities (best split):", max(membership(x)), "\n")
    cat("Modularity (best split):", modularity(x), "\n")
    cat("Membership vector:\n")
    print(membership(x))
  } else if (algorithm(x) == "optimal") {
    cat("Number of communities:", max(membership(x)), "\n")
    cat("Modularity:", modularity(x), "\n")
    cat("Membership vector:\n")
    print(membership(x))
  } else if (algorithm(x) == "infomap") {
    cat("Number of communities:", max(membership(x)), "\n")
    if (!is.null(x$modularity)) {
      cat("Modularity:", modularity(x), "\n")
    }
    cat("Membership vector:\n")
    print(membership(x))
  } else {
    cat("Number of communities:", max(membership(x)), "\n")
    if (!is.null(x$modularity)) {
      cat("Modularity:", modularity(x), "\n")
    }
    cat("Membership vector:\n")
    print(membership(x))
  }
  invisible(x)
}

create.communities <- function(membership, algorithm=NULL, merges=NULL,
                               modularity=NULL, ...) {

  stopifnot(is.numeric(membership))
  stopifnot(is.null(algorithm) ||
            (is.character(algorithm) && length(algorithm)==1))
  stopifnot(is.null(merges) ||
            (is.matrix(merges) && is.numeric(merges) && ncol(merges)==2))
  stopifnot(is.null(modularity) ||
            (is.numeric(modularity) &&
             length(modularity) %in% c(1, length(membership))))

  res <- list(membership=membership,
              algorithm=if (is.null(algorithm)) "unknown" else algorithm,
              modularity=modularity, ...)
  if (!is.null(merges)) {
    res$merges <- merges
  }
  class(res) <- "communities"
  res
}

modularity <- function(x, ...)
  UseMethod("modularity")

modularity.igraph <- function(x, membership, weights=NULL, ...) {
  # Argument checks
  if (!is.igraph(x)) { stop("Not a graph object") }
  membership <- as.numeric(membership)
  if (!is.null(weights)) weights <- as.numeric(weights)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_modularity", x, membership-1, weights,
        PACKAGE="igraph")
  res
}

modularity.communities <- function(x, ...) {
  if (!is.null(x$modularity)) {
    max(x$modularity)
  } else {
    stop("Modularity was not calculated")
  }
}

length.communities <- function(x) {
  m <- membership(x)
  max(m)
}

sizes <- function(communities) {
  m <- membership(communities)
  table(`Community sizes`=m)
}

communities <- function(communities) {
  m <- membership(communities)
  tapply(seq_along(m), m, simplify=FALSE,
         function(x) x)
}

algorithm <- function(communities) {
  communities$algorithm
}

merges <- function(communities) {
  if (!is.null(communities$merges)) {
    communities$merges
  } else {
    stop("Not a hierarchical community structure")
  }
}

crossing <- function(communities, graph) {
  m <- membership(communities)
  el <- get.edgelist(graph, names=FALSE)
  m1 <- m[el[,1]]
  m2 <- m[el[,2]]
  res <- m1 != m2
  if (!is.null(names(m1))) {
    names(res) <- paste(names(m1), names(m2), sep="|")
  }
  res
}

code.length <- function(communities) {
  communities$codelength
}

is.hierarchical <- function(communities, full=FALSE) {
  alg <- algorithm(communities)
  if (alg %in% c("walktrap", "edge betweenness","fast greedy") ||
      (alg == "leading eigenvector" && !full)) {
    TRUE
  } else if (alg %in% c("spinglass", "label propagation", "multi level",
                        "optimal") ||
             (alg == "leading eigenvector" && full)) {
    FALSE
  } else {
    stop("Unknown community detection algorithm")
  }
}

complete.dend <- function(comm, use.modularity) {
  merges <- comm$merges
  if (nrow(merges) < comm$vcount-1) {
    if (use.modularity) {
      stop(paste("`use.modularity' requires a full dendrogram,",
                 "i.e. a connected graph"))
    }
    miss <- seq_len(comm$vcount + nrow(merges))[-as.vector(merges)]
    miss <- c(miss, seq_len(length(miss)-2) + comm$vcount+nrow(merges))
    miss <- matrix(miss, byrow=TRUE, ncol=2)
    merges <- rbind(merges, miss)
  }
  storage.mode(merges) <- "integer"

  merges
}

# The following functions were adapted from the stats R package

as.dendrogram.communities <- function(object, hang=-1, use.modularity=FALSE,
                                      ...) {
  if (!is.hierarchical(object, full=TRUE)) {
    stop("Not a fully hierarchical community structure")
  }

  .memberDend <- function(x) {
    r <- attr(x,"x.member")
    if(is.null(r)) {
      r <- attr(x,"members")
    if(is.null(r)) r <- 1:1
    }
    r
  }

  ## If multiple components, then we merge them in arbitrary order
  merges <- complete.dend(object, use.modularity)
  
  storage.mode(merges) <- "integer"
  
  if (is.null(object$names)) {
    object$names <- 1:(nrow(merges)+1)
  }
  z <- list()
  if (!use.modularity || is.null(object$modularity)) {
    object$height <- 1:nrow(merges)
  } else {
    object$height <- object$modularity[-1]
    object$height <- cumsum(object$height - min(object$height))
  }
  nMerge <- length(oHgt <- object$height)
  if (nMerge != nrow(merges))
    stop("'merge' and 'height' do not fit!")
  hMax <- oHgt[nMerge]
  one <- 1L
  two <- 2L
  leafs <- nrow(merges)+1
  for (k in 1:nMerge) {
    x <- merges[k, ]# no sort() anymore!
    if (any(neg <- x < leafs+1))
      h0 <- if (hang < 0) 0 else max(0, oHgt[k] - hang * hMax)
    if (all(neg)) {                  # two leaves
      zk <- as.list(x)
      attr(zk, "members") <- two
      attr(zk, "midpoint") <- 0.5 # mean( c(0,1) )
      objlabels <- object$names[x]
      attr(zk[[1]], "label") <- objlabels[1]
      attr(zk[[2]], "label") <- objlabels[2]
      attr(zk[[1]], "members") <- attr(zk[[2]], "members") <- one
      attr(zk[[1]], "height") <- attr(zk[[2]], "height") <- h0
      attr(zk[[1]], "leaf") <- attr(zk[[2]], "leaf") <- TRUE
    }
    else if (any(neg)) {            # one leaf, one node
      X <- as.character(x)
      ## Originally had "x <- sort(..) above => leaf always left, x[1];
      ## don't want to assume this
      isL <- x[1] < leafs+1 ## is leaf left?
      zk <-
        if(isL) list(x[1], z[[X[2]]])
        else    list(z[[X[1]]], x[2])
      attr(zk, "members") <- attr(z[[X[1 + isL]]], "members") + one
      attr(zk, "midpoint") <-
        (.memberDend(zk[[1]]) + attr(z[[X[1 + isL]]], "midpoint"))/2
      attr(zk[[2 - isL]], "members") <- one
      attr(zk[[2 - isL]], "height") <- h0
      attr(zk[[2 - isL]], "label") <- object$names[x[2 - isL]]
      attr(zk[[2 - isL]], "leaf") <- TRUE
      }
    else {                        # two nodes
      x <- as.character(x)
      zk <- list(z[[x[1]]], z[[x[2]]])
      attr(zk, "members") <- attr(z[[x[1]]], "members") +
        attr(z[[x[2]]], "members")
      attr(zk, "midpoint") <- (attr(z[[x[1]]], "members") +
                               attr(z[[x[1]]], "midpoint") +
                               attr(z[[x[2]]], "midpoint"))/2
    }
    attr(zk, "height") <- oHgt[k]
    z[[k <- as.character(k+leafs)]] <- zk
  }
  z <- z[[k]]
  class(z) <- "dendrogram"
  z
}

as.hclust.communities <- function(x, hang=-1, use.modularity=FALSE,
                                  ...) {
  as.hclust(as.dendrogram(x, hang=hang, use.modularity=use.modularity))
}

asPhylo <- function(x, ...)
  UseMethod("asPhylo")

asPhylo.communities <- function(x, use.modularity=FALSE, ...) {

  if (!is.hierarchical(x, full=TRUE)) {
    stop("Not a fully hierarchical community structure")
  }

  require(ape, quietly = TRUE)
  
  ## If multiple components, then we merge them in arbitrary order
  merges <- complete.dend(x, use.modularity)

  if (!use.modularity || is.null(x$modularity)) {
    height <- 1:nrow(merges)
  } else {
    height <- x$modularity[-1]
    height <- cumsum(height - min(height))
  }

  if (is.null(x$names)) {
    labels <- 1:(nrow(merges)+1)
  } else {
    labels <- x$names
  }

  N <- nrow(merges)
  edge <- matrix(0L, 2*N, 2)
  edge.length <- numeric(2*N)
  node <- integer(N)
  node[N] <- N + 2L
  cur.nod <- N + 3L
  j <- 1L
  for (i in N:1) {
    edge[j:(j+1), 1] <- node[i]
    for (l in 1:2) {
      k <- j + l -1L
      y <- merges[i, l]
      if (y > N+1) {
        edge[k, 2] <- node[y-N-1] <- cur.nod
        cur.nod <- cur.nod + 1L
        edge.length[k] <- height[i] - height[y-N-1]
      } else {
        edge[k, 2] <- y
        edge.length[k] <- height[i]
      }
    }
    j <- j + 2L    
  }

  obj <- list(edge=edge, edge.length=edge.length/2, tip.label=labels,
              Nnode=N)
  class(obj) <- "phylo"
  reorder(obj)
}

cutat <- function(communities, no, steps) {

  if (!inherits(communities, "communities")) {
    stop("Not a community structure")
  }
  if (!is.hierarchical(communities, full=TRUE)) {
    stop("Not a fully hierarchical communitity structure")
  }

  if ((!missing(no) && !missing(steps)) ||
      ( missing(no) &&  missing(steps))) {
    stop("Please give either `no' or `steps' (but not both)")
  }

  if (!missing(steps)) {
    mm <- merges(communities)
    if (steps > nrow(mm)) {
      warning("Cannot make that many steps")
      steps <- nrow(mm)
    }
    community.to.membership2(mm, communities$vcount, steps)
  } else {
    mm <- merges(communities)
    noc <- communities$vcount - nrow(mm) # final number of communities
    if (no<noc) {
      warning("Cannot have that few communities")
      no=noc
    }
    steps <- communities$vcount-no
    community.to.membership2(mm, communities$vcount, steps)    
  }
}

showtrace <- function(communities) {

  if (!inherits(communities, "communities")) {
    stop("Not a community structure")
  }
  if (is.null(communities$history)) {
    stop("History was not recorded")
  }

  res <- character()
  i <- 1
  while (i <= length(communities$history)) {
    if (communities$history[i] == 2) {  # IGRAPH_LEVC_HIST_SPLIT
      resnew <- paste("Splitting community", communities$history[i+1],
                      "into two.")
      i <- i + 2
    } else if (communities$history[i]==3) { # IGRAPH_LEVC_HIST_FAILED
      resnew <- paste("Failed splitting community",
                      communities$history[i+1], "into two.")
      i <- i + 2
    } else if (communities$history[i]==4) { # IGRAPH_LEVC_START_FULL
      resnew <- "Starting with the whole graph as a community."
      i <- i + 1
    } else if (communities$history[i]==5) { # IGRAPH_LEVC_START_GIVEN
      resnew <- paste("Starting from the", communities$history[i+1],
                      "given communities.")
      i <- i + 2
    }

    res <- c(res, resnew)
  }
  res
}

#####################################################################

community.to.membership2 <- function(merges, vcount, steps) {
  mode(merges) <- "numeric"
  mode(vcount) <- "numeric"
  mode(steps)  <- "numeric"
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_community_to_membership2", merges-1, vcount, steps,
               PACKAGE="igraph")
  res+1
}

#####################################################################

spinglass.community <- function(graph, weights=NULL, vertex=NULL, spins=25,
                                parupdate=FALSE, start.temp=1,
                                stop.temp=0.01, cool.fact=0.99,
                                update.rule=c("config", "random", "simple"),
                                gamma=1.0, implementation=c("orig", "neg"),
                                gamma.minus=1.0) {

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

  update.rule <- igraph.match.arg(update.rule)
  update.rule <- switch(update.rule, "simple"=0, "random"=0, "config"=1)
  implementation <- switch(igraph.match.arg(implementation),
                                            "orig"=0, "neg"=1)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  if (is.null(vertex)) {    
    res <- .Call("R_igraph_spinglass_community", graph, weights,
                 as.numeric(spins), as.logical(parupdate),
                 as.numeric(start.temp),
                 as.numeric(stop.temp), as.numeric(cool.fact),
                 as.numeric(update.rule), as.numeric(gamma),
                 as.numeric(implementation), as.numeric(gamma.minus),
                 PACKAGE="igraph")
    res$algorithm  <- "spinglass"
    res$vcount     <- vcount(graph)
    res$membership <- res$membership + 1
    if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
      res$names <- get.vertex.attribute(graph, "name")
    }
    class(res) <- "communities"
  } else {
    res <- .Call("R_igraph_spinglass_my_community", graph, weights,
                 as.igraph.vs(graph, vertex)-1, as.numeric(spins), 
                 as.numeric(update.rule), as.numeric(gamma),
                 PACKAGE="igraph")
    res$community <- res$community + 1
  }
  res
}

walktrap.community <- function(graph, weights=E(graph)$weight, steps=4,
                               merges=TRUE, modularity=TRUE,
                               membership=TRUE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  if (membership && !modularity) {
    modularity <- TRUE
  }
  
  if (!is.null(weights)) {
    weights <- as.numeric(weights)
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_walktrap_community", graph, weights, as.numeric(steps),
        as.logical(merges), as.logical(modularity), as.logical(membership),
        PACKAGE="igraph")
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    res$names <- V(graph)$name
  }

  res$vcount <- vcount(graph)
  res$algorithm <- "walktrap"
  res$membership <- res$membership + 1
  res$merges <- res$merges + 1
  class(res) <- "communities"
  res
}

edge.betweenness.community <- function(graph, weights=E(graph)$weight,
                                       directed=TRUE,
                                       edge.betweenness=TRUE,
                                       merges=TRUE, bridges=TRUE,
                                       modularity=TRUE,
                                       membership=TRUE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  if (!is.null(weights)) {
    weights <- as.numeric(weights)
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_community_edge_betweenness", graph, weights,
               as.logical(directed),
               as.logical(edge.betweenness),
               as.logical(merges), as.logical(bridges),
               as.logical(modularity), as.logical(membership),
               PACKAGE="igraph")
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    res$names <- V(graph)$name
  }
  res$vcount <- vcount(graph)
  res$algorithm <- "edge betweenness"
  res$membership <- res$membership + 1
  res$merges <- res$merges + 1
  res$removed.edges <- res$removed.edges + 1
  res$bridges <- res$bridges + 1
  class(res) <- "communities"
  res
}

edge.betweenness.community.merges <- function(graph, edges) {
  
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_community_eb_get_merges", graph, as.numeric(edges),
        PACKAGE="igraph")
}

fastgreedy.community <- function(graph, merges=TRUE, modularity=TRUE,
                                 membership=TRUE,
                                 weights=E(graph)$weight) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (!is.null(weights)) {
    weights <- as.numeric(weights)
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_community_fastgreedy", graph, as.logical(merges),
               as.logical(modularity), as.logical(membership), weights,
               PACKAGE="igraph")
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    res$names <- V(graph)$name
  }
  res$algorithm <- "fast greedy"
  res$vcount <- vcount(graph)
  res$membership <- res$membership + 1
  res$merges <- res$merges + 1
  class(res) <- "communities"
  res
}

community.to.membership <- function(graph, merges, steps, membership=TRUE,
                                    csize=TRUE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  merges <- as.matrix(merges)
  merges <- structure(as.numeric(merges), dim=dim(merges))
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_community_to_membership", graph, merges-1,
        as.numeric(steps), as.logical(membership), as.logical(csize),
        PACKAGE="igraph")
}

igraph.i.levc.arp <- function(externalP, externalE) {
  f <- function(v) {
    v <- as.numeric(v)
    base::.Call("R_igraph_i_levc_arp", externalP, externalE, v,
                PACKAGE="igraph");
  }
  f
}

leading.eigenvector.community <- function(graph, steps=-1, weights=NULL,
                                          start=NULL,
                                          options=igraph.arpack.default,
                                          callback=NULL, extra=NULL,
                                          env=parent.frame()){

  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  steps <- as.integer(steps)
  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) { 
    weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
    weights <- as.numeric(weights) 
  } else { 
    weights <- NULL 
  }
  if (!is.null(start)) { start <- as.numeric(start)-1 }
  options.tmp <- igraph.arpack.default; options.tmp[ names(options) ] <- options ; options <- options.tmp
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_community_leading_eigenvector", graph, steps,
               weights, options, start, callback, extra, env,
               environment(igraph.i.levc.arp),
               PACKAGE="igraph")
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    res$names <- V(graph)$name
  }
  res$algorithm <- "leading eigenvector"
  res$vcount <- vcount(graph)
  res$membership <- res$membership + 1
  res$merges <- res$merges + 1
  res$history <- res$history + 1
  class(res) <- "communities"
  res
}

label.propagation.community <- function(graph, weights=NULL, initial=NULL,
                                        fixed=NULL) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) { 
  weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
  weights <- as.numeric(weights) 
  } else { 
  weights <- NULL 
  }
  if (!is.null(initial)) initial <- as.numeric(initial)
  if (!is.null(fixed)) fixed <- as.logical(fixed)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_community_label_propagation", graph, weights, initial, fixed,
        PACKAGE="igraph")
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    res$names <- V(graph)$name
  }
  res$vcount <- vcount(graph)
  res$algorithm <- "label propagation"
  res$membership <- res$membership + 1
  class(res) <- "communities"
  res
}

multilevel.community <- function(graph, weights=NULL) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
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
  res <- .Call("R_igraph_community_multilevel", graph, weights,
        PACKAGE="igraph")
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    res$names <- V(graph)$name
  }
  res$vcount <- vcount(graph)
  res$algorithm <- "multi level"
  res$membership <- res$membership + 1
  res$memberships <- res$memberships + 1
  class(res) <- "communities"
  res
}

optimal.community <- function(graph, weights=NULL) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
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
  res <- .Call("R_igraph_community_optimal_modularity", graph, weights,
               PACKAGE="igraph")
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    res$names <- V(graph)$name
  }
  res$vcount <- vcount(graph)
  res$algorithm <- "optimal"
  res$membership <- res$membership + 1
  class(res) <- "communities"
  res
}

infomap.community <- function(graph, e.weights=NULL, v.weights=NULL,
                              nb.trials=10, modularity=TRUE) {
  
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  if (is.null(e.weights) && "weight" %in% list.edge.attributes(graph)) { 
    e.weights <- E(graph)$weight 
  } 
  if (!is.null(e.weights) && any(!is.na(e.weights))) { 
    e.weights <- as.numeric(e.weights) 
  } else { 
    e.weights <- NULL 
  }
  if (is.null(v.weights) && "weight" %in% list.vertex.attributes(graph)) { 
    v.weights <- V(graph)$weight 
  } 
  if (!is.null(v.weights) && any(!is.na(v.weights))) { 
    v.weights <- as.numeric(v.weights) 
  } else { 
    v.weights <- NULL 
  }
  nb.trials <- as.integer(nb.trials)
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_community_infomap", graph, e.weights,
               v.weights, nb.trials,
               PACKAGE="igraph")

  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    res$names <- V(graph)$name
  }
  res$vcount <- vcount(graph)
  res$algorithm <- "infomap"
  res$membership <- res$membership + 1
  if (modularity) {
    res$modularity <- modularity(graph, res$membership, weights=e.weights)
  }
  class(res) <- "communities"
  res
}

plot.communities <- function(x, y,
                             colbar=rainbow(length(x)),
                             col=colbar[membership(x)],
                             mark.groups=communities(x),
                             edge.color=c("black", "red")[crossing(x,y)+1],
                             ...) {

  plot(y, vertex.color=col, mark.groups=mark.groups,
       edge.color=edge.color,
       ...)  
}

dendPlot <- function(x, mode=getIgraphOpt("dend.plot.type"), ...)
  UseMethod("dendPlot")

dendPlot.communities <- function(x, 
                                 mode=getIgraphOpt("dend.plot.type"), ...,
                                 use.modularity=FALSE) {  
  mode <- igraph.match.arg(mode, c("auto", "phylo", "hclust", "dendrogram"))

  if (mode=="auto") {
    value <- tryCatch(suppressWarnings(library("ape", character.only=TRUE,
                                               logical.return=TRUE,
                                               warn.conflicts=FALSE,
                                               quietly=TRUE,
                                               pos="package:base")),
                      error=function(e) e)
    mode <- if (value) "phylo" else "hclust"
  }
  
  if (mode=="hclust") {
    dendPlotHclust(x, use.modularity=use.modularity, ...)
  } else if (mode=="dendrogram") {
    dendPlotDendrogram(x, use.modularity=use.modularity, ...)
  } else if (mode=="phylo") {
    dendPlotPhylo(x, use.modularity=use.modularity, ...)
  }
}

dendPlotHclust <- function(communities, rect=length(communities),
                           colbar=rainbow(rect), hang=-1, ann=FALSE,
                           main="", sub="", xlab="", ylab="", ...,
                           use.modularity=FALSE) {
  hc <- as.hclust(communities, hang=hang, use.modularity=use.modularity)
  ret <- plot(hc, hang=hang, ann=ann, main=main, sub=sub, xlab=xlab,
              ylab=ylab, ...)
  if (rect > 0) {
    rect.hclust(hc, k=rect, border=colbar)
  }
  invisible(ret)
}

dendPlotDendrogram <- function(communities, hang=-1, ...,
                               use.modularity=FALSE) {
  plot(as.dendrogram(communities, hang=hang, use.modularity=use.modularity),
       ...)
}

dendPlotPhylo <- function(communities, colbar=rainbow(length(communities)),
                          col=colbar[membership(communities)],
                          mark.groups=communities(communities),
                          use.modularity=FALSE, 
                          edge.color="#AAAAAAFF",
                          edge.lty=c(1,2), ...) {
  
  phy <- asPhylo(communities, use.modularity=use.modularity)

  getedges <- function(tip) {
    repeat {      
      ee <- which(! phy$edge[,1] %in% tip & phy$edge[,2] %in% tip)
      if (length(ee)<=1) { break }
      tip <- c(tip, unique(phy$edge[ee,1]))
    }
    ed <- which(phy$edge[,1] %in% tip & phy$edge[,2] %in% tip)
    eds <- phy$edge[ed, 1]
    good <- which(phy$edge[ed,1] %in% which(tabulate(eds) != 1))
    ed[good]
  }
  gredges <- lapply(mark.groups, getedges)

  if (length(mark.groups) > 0) {
    ecol <- rep(edge.color, nrow(phy$edge))
    for (gr in seq_along(gredges)) {
      ecol[gredges[[gr]]] <- colbar[gr]
    }
  } else {
    ecol <- edge.color
  }
  
  elty <- rep(edge.lty[2], nrow(phy$edge))
  elty[ unlist(gredges) ] <- edge.lty[1]
  
  plot(phy, edge.color=ecol, edge.lty=elty, tip.color=col, ...)
}

compare <- function(comm1, comm2, method=c("vi", "nmi",
                                       "split.join", "rand",
                                       "adjusted.rand"))
  UseMethod("compare")

compare.communities <- function(comm1, comm2, method=c("vi", "nmi",
                                                "split.join", "rand",
                                                "adjusted.rand")) {
  compare.numeric(comm1, comm2, method)
}

compare.numeric <- function(comm1, comm2, method=c("vi", "nmi",
                                            "split.join", "rand",
                                            "adjusted.rand")) {
  comm1 <- if (inherits(comm1, "communities")) {
    membership(comm1)
  } else {
    as.numeric(comm1)
  }
  comm2 <- if (inherits(comm2, "communities")) {
    membership(comm2)
  } else {
    as.numeric(comm2)
  }
  method <- switch(igraph.match.arg(method), vi = 0, nmi = 1, 
                   split.join = 2, rand = 3, adjusted.rand = 4)
  on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
  res <- .Call("R_igraph_compare_communities", comm1, comm2, 
               method, PACKAGE = "igraph")
  res  
}
                                
compare.default <- function(comm1, comm2, method=c("vi", "nmi",
                                            "split.join", "rand",
                                            "adjusted.rand")) {
  compare.numeric(as.numeric(comm1), as.numeric(comm2), method)
}

