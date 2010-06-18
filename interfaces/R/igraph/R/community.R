
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

membership <- function(commmunities, ...)
  UseMethod("membership")

membership.communities <- function(communities, ...) {
  if (!is.null(communities$membership)) {
    return(communities$membership)
  } else if (!is.null(communities$merges) &&
             !is.null(communities$modularity)) {
    return(community.to.membership2(communities$merges, communities$vcount,
                                    which.max(communities$modularity)-1))
  } else {
    stop("Cannot calculate community membership")
  }
}

print.communities <- function(x, ...) {
  cat("Graph community structure calculated with the",
      algorithm(x), "algorithm\n")
  if (algorithm(x)=="spinglass") {
    cat("Number of communities:", max(membership(x))+1, "\n")
    cat("Modularity:", modularity(x), "\n")
    cat("Membership vector:\n")
    print(membership(x))
  } else if (algorithm(x) %in% c("walktrap", "edge betweenness",
                                "fast greedy")) {
    cat("Number of communities (best split):", max(membership(x))+1, "\n")
    cat("Modularity (best split):", modularity(x), "\n")
    cat("Membership vector:\n")
    print(membership(x))
  } else if (algorithm(x) %in% c("leading eigenvector",
                                "leading eigenvector, naive")) {
    cat("Number of communities (best split):", max(membership(x))+1, "\n")
    cat("Modularity (best split):", modularity(x), "\n")
    cat("Membership vector:\n")
    print(membership(x))
  } else if (algorithm(x) == "label propagation") {
    cat("Number of communities:", max(membership(x))+1, "\n")
    cat("Modularity:", modularity(x), "\n")
    cat("Membership vector:\n")
    print(membership(x))
  } else if (algorithm(x) == "multi level") {
    cat("Number of communities (best split):", max(membership(x))+1, "\n")
    cat("Modularity (best split):", modularity(x), "\n")
    cat("Membership vector:\n")
    print(membership(x))
  } else if (algorithm(x) == "optimal") {
    cat("Number of communities:", max(membership(x))+1, "\n")
    cat("Modularity:", modularity(x), "\n")
    cat("Membership vector:\n")
    print(membership(x))
  }    
}

modularity <- function(x, ...)
  UseMethod("modularity")

modularity.communities <- function(communities, ...) {
  if (!is.null(communities$modularity)) {
    max(communities$modularity)
  } else {
    stop("Modularity was not calculated")
  }
}

length.communities <- function(x) {
  m <- membership(x)
  max(m)+1
}

sizes <- function(x, ...)
  UseMethod("sizes")

sizes.communities <- function(x, ...) {
  m <- membership(x)
  table(`Community sizes`=m)
}

algorithm <- function(x, ...)
  UseMethod("algorithm")

algorithm.communities <- function(x, ...) {
  x$algorithm
}

merges <- function(x, ...)
  UseMethod("merges")

merges.communities <- function(x, ...) {
  if (!is.null(x$merges)) {
    x$merges
  } else {
    stop("Not a hierarchical community structure")
  }
}

is.hierarchical <- function(x, ...)
  UseMethod("is.hierarchical")

is.hierarchical.communities <- function(x, ...) {
  if (algorithm(x) %in% c("walktrap", "edge betweenness", "fast greedy")) {
    TRUE
  } else if (algorithm(x) %in% c("spinglass", "leading eigenvector",
                                 "leading eigenvector, naive",
                                 "label propagation", "multi level",
                                 "optimal")) {
    FALSE
  } else {
    stop("Unknown community detection algorithm")
  }
}

# The following functions were adapted from the stats R package

as.dendrogram.communities <- function(object, hang=-1,
                                      use.modularity=FALSE, ...) {
  if (!is.hierarchical(object)) {
    stop("Not a hierarchical community structure")
  }

  .memberDend <- function(x) {
    r <- attr(x,"x.member")
    if(is.null(r)) {
      r <- attr(x,"members")
    if(is.null(r)) r <- 1:1
    }
    r
  }
  
  if (is.null(object$labels))
    object$labels <- 1:(nrow(object$merges)+1)-1
  z <- list()
  if (!use.modularity || is.null(object$modularity)) {
    object$height <- 1:nrow(object$merges)
  } else {
    object$height <- object$modularity[-1]
    object$height <- cumsum(object$height - min(object$height))
  }
  nMerge <- length(oHgt <- object$height)
  if (nMerge != nrow(object$merges))
    stop("'merge' and 'height' do not fit!")
  hMax <- oHgt[nMerge]
  one <- 1:1;
  two <- 2:2 # integer!
  leafs <- nrow(object$merges)+1
  for (k in 1:nMerge) {
    x <- object$merges[k, ]# no sort() anymore!
    if (any(neg <- x < leafs))
      h0 <- if (hang < 0) 0 else max(0, oHgt[k] - hang * hMax)
    if (all(neg)) {                  # two leaves
      zk <- as.list(x)
      attr(zk, "members") <- two
      attr(zk, "midpoint") <- 0.5 # mean( c(0,1) )
      objlabels <- object$labels[x+1]
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
      isL <- x[1] < leafs ## is leaf left?
      zk <-
        if(isL) list(x[1], z[[X[2]]])
        else    list(z[[X[1]]], x[2])
      attr(zk, "members") <- attr(z[[X[1 + isL]]], "members") + one
      attr(zk, "midpoint") <-
        (.memberDend(zk[[1]]) + attr(z[[X[1 + isL]]], "midpoint"))/2
      attr(zk[[2 - isL]], "members") <- one
      attr(zk[[2 - isL]], "height") <- h0
      attr(zk[[2 - isL]], "label") <- object$labels[x[2 - isL]+1]
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
    z[[k <- as.character(k+leafs-1)]] <- zk
  }
  z <- z[[k]]
  class(z) <- "dendrogram"
  z
}

cutat <- function(communities, no, steps) {

  if (!inherits(communities, "communities")) {
    stop("Not a community structure")
  }
  if (!is.hierarchical(communities)) {
    stop("Not a hierarchical communitity structure")
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

#####################################################################

community.to.membership2 <- function(merges, vcount, steps) {
  mode(merges) <- "numeric"
  mode(vcount) <- "numeric"
  mode(steps)  <- "numeric"
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_community_to_membership2", merges, vcount, steps,
        PACKAGE="igraph")
}

#####################################################################

spinglass.community <- function(graph, weights=NULL, vertex=NULL, spins=25,
                                parupdate=FALSE, start.temp=1,
                                stop.temp=0.01, cool.fact=0.99,
                                update.rule=c("config", "random", "simple"),
                                gamma=1.0, implementation=c("orig", "neg"),
                                lambda=1.0) {

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
                 as.numeric(implementation), as.numeric(lambda),
                 PACKAGE="igraph")
    res$algorithm <- "spinglass"
    res$vcount    <- vcount(graph)
    class(res) <- "communities"
  } else {
    res <- .Call("R_igraph_spinglass_my_community", graph, weights,
                 as.igraph.vs(graph, vertex), as.numeric(spins), 
                 as.numeric(update.rule), as.numeric(gamma),
                 PACKAGE="igraph")
  }
  res
}

walktrap.community <- function(graph, weights=E(graph)$weight, steps=4,
                               merges=TRUE, modularity=TRUE, labels=TRUE,
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
  if (labels && "name" %in% list.vertex.attributes(graph)) {
    res$labels <- V(graph)$name
  }

  res$vcount <- vcount(graph)
  res$algorithm <- "walktrap"
  class(res) <- "communities"
  res
}

edge.betweenness.community <- function(graph, directed=TRUE,
                                       edge.betweenness=TRUE,
                                       merges=TRUE, bridges=TRUE,
                                       labels=TRUE, modularity=TRUE,
                                       membership=TRUE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_community_edge_betweenness", graph, as.logical(directed),
               as.logical(edge.betweenness),
               as.logical(merges), as.logical(bridges),
               as.logical(modularity), as.logical(membership),
               PACKAGE="igraph")
  if (labels && "name" %in% list.vertex.attributes(graph)) {
    res$labels <- V(graph)$name
  }
  res$vcount <- vcount(graph)
  res$algorithm <- "edge betweenness"
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
  res$algorithm <- "fast greedy"
  res$vcount <- vcount(graph)
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
  .Call("R_igraph_community_to_membership", graph, merges, as.numeric(steps),
        as.logical(membership), as.logical(csize),
        PACKAGE="igraph")
}

leading.eigenvector.community <- function(graph, steps=-1, options=igraph.arpack.default) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  steps <- as.numeric(steps)
  options.tmp <- igraph.arpack.default; options.tmp[ names(options) ] <- options ; options <- options.tmp

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_community_leading_eigenvector", graph, steps, options,
        PACKAGE="igraph")
  res$algorithm <- "leading eigenvector"
  res$vcount <- vcount(graph)
  class(res) <- "communities"
  res
}

leading.eigenvector.community.naive <- function(graph, steps=-1, options=igraph.arpack.default) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  steps <- as.numeric(steps)
  options.tmp <- igraph.arpack.default; options.tmp[ names(options) ] <- options ; options <- options.tmp

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_community_leading_eigenvector_naive", graph, steps, options,
        PACKAGE="igraph")
  res$algorithm <- "leading eigenvector, naive"
  res$vcount <- vcount(graph)
  class(res) <- "communities"
  res
}

leading.eigenvector.community.step <- function(graph, fromhere=NULL,
                                               membership=rep(0, vcount(graph)),
                                               community=0,
                                               options=igraph.arpack.default) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }
  if (!is.null(fromhere)) {
    if (class(fromhere) != "igraph.eigencstep") {
      stop("invalid community structure object given")
    }
    membership <- fromhere[["membership"]]
  }

  options.tmp <- igraph.arpack.default
  options.tmp[names(options)] <- options
  options <- options.tmp

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_community_leading_eigenvector_step",
               graph, as.numeric(membership), as.numeric(community),
               options,
               PACKAGE="igraph")
  class(res) <- "igraph.eigencstep"
  res
}

label.propagation.community <- function(graph, weights=NULL, initial=NULL, fixed=NULL) {
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
  res$vcount <- vcount(graph)
  res$algorithm <- "label propagation"
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
  res$vcount <- vcount(graph)
  res$algorithm <- "multi level"
  class(res) <- "communities"
  res
}

optimal.community <- function(graph, verbose=igraph.par("verbose")) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  verbose <- as.logical(verbose)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_community_optimal_modularity", graph, verbose,
        PACKAGE="igraph")
  res$vcount <- vcount(graph)
  res$algorithm <- "optimal"
  class(res) <- "communities"
  res
}

plot.communities <- function(communities, graph,
                             colbar=rainbow(length(communities), ...) {

  col <- colbar[membership(communities)+1]
  plot(graph, vertex.color=col, ...)  
}
