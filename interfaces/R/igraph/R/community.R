
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
                                update.rule=c("config", "random", "simple"),
                                gamma=1.0) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(weights)) {
    if ("weight" %in% list.edge.attributes(graph)) {
      weights <- as.numeric(E(graph)$weight)
    } else {
      weights <- as.numeric(rep(1, ecount(graph)))
    }
  }

  update.rule <- igraph.match.arg(update.rule)
  update.rule <- switch(update.rule, "simple"=0, "random"=0, "config"=1)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  if (is.null(vertex)) {    
    .Call("R_igraph_spinglass_community", graph, weights,
          as.numeric(spins), as.logical(parupdate), as.numeric(start.temp),
          as.numeric(stop.temp), as.numeric(cool.fact),
          as.numeric(update.rule), as.numeric(gamma),
          PACKAGE="igraph0")
  } else {
    .Call("R_igraph_spinglass_my_community", graph, weights,
          as.igraph.vs(graph, vertex), as.numeric(spins), 
          as.numeric(update.rule), as.numeric(gamma),
          PACKAGE="igraph0")
  }    
}

walktrap.community <- function(graph, weights=E(graph)$weight, steps=4, merges=TRUE,
                               modularity=TRUE, labels=TRUE, membership=TRUE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  if (membership && !modularity) {
    modularity <- TRUE
  }
  
  if (!is.null(weights)) {
    weights <- as.numeric(weights)
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  res <- .Call("R_igraph_walktrap_community", graph, weights, as.numeric(steps),
        as.logical(merges), as.logical(modularity),
        PACKAGE="igraph0")
  if (labels && "name" %in% list.vertex.attributes(graph)) {
    res$labels <- V(graph)$name
  }

  res <- append(res, list(membership=NULL))
  if (membership) {
    res$membership <-
      community.to.membership(graph, res$merges,
                              steps=which.max(res$modularity)-1)$membership
  }
  
  class(res) <- "igraph.walktrap"
  res
}

# The following two functions were adapted from the stats R package

.memberDend <- function(x) {
  r <- attr(x,"x.member")
  if(is.null(r)) {
    r <- attr(x,"members")
    if(is.null(r)) r <- 1:1
  }
  r
}

as.dendrogram.igraph.walktrap <- function (object, hang=-1,
                                           use.modularity=FALSE, ...)
{
  stopifnot(nrow(object$merges)> 0)
  storage.mode(object$merges) <- "integer"
  object$merges <- object$merges + 1L
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
    if (any(neg <- x < leafs+1))
      h0 <- if (hang < 0) 0 else max(0, oHgt[k] - hang * hMax)
    if (all(neg)) {                  # two leaves
      zk <- as.list(x)
      attr(zk, "members") <- two
      attr(zk, "midpoint") <- 0.5 # mean( c(0,1) )
      objlabels <- object$labels[x]
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
      attr(zk[[2 - isL]], "label") <- object$labels[x[2 - isL]]
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

edge.betweenness.community <- function(graph, directed=TRUE,
                                       edge.betweenness=TRUE,
                                       merges=TRUE, bridges=TRUE,
                                       labels=TRUE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  res <- .Call("R_igraph_community_edge_betweenness", graph, as.logical(directed),
               as.logical(edge.betweenness),
               as.logical(merges), as.logical(bridges),
               PACKAGE="igraph0")
  if (labels && "name" %in% list.vertex.attributes(graph)) {
    res$labels <- V(graph)$name
  }
  class(res) <- "igraph.ebc"
  res
}

edge.betweenness.community.merges <- function(graph, edges) {
  
  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_community_eb_get_merges", graph, as.numeric(edges),
        PACKAGE="igraph0")
}

# Adapted from the stats package

as.dendrogram.igraph.ebc <- function (object, hang=-1,
                                      use.eb=FALSE,
                                      merges=NULL,
                                      bridges=NULL,
                                      ...)
{
  if (!is.null(merges)) {
    object$merges <- merges
  }
  stopifnot(nrow(object$merges)> 0)
  storage.mode(object$merges) <- "integer"
  object$merges <- object$merges + 1L
  if (is.null(object$labels))
    object$labels <- 1:(nrow(object$merges)+1)-1
  z <- list()
  if (!use.eb || is.null(object$edge.betweenness) || is.null(object$bridges)) {
    object$height <- 1:nrow(object$merges)
  } else {
    object$height <- object$edge.betweenness[object$bridges]
    object$height <- cumsum(object$height)
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
    if (any(neg <- x < leafs+1))
      h0 <- if (hang < 0) 0 else max(0, oHgt[k] - hang * hMax)
    if (all(neg)) {                  # two leaves
      zk <- as.list(x)
      attr(zk, "members") <- two
      attr(zk, "midpoint") <- 0.5 # mean( c(0,1) )
      objlabels <- object$labels[x]
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
      attr(zk[[2 - isL]], "label") <- object$labels[x[2 - isL]]
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

fastgreedy.community <- function(graph, merges=TRUE, modularity=TRUE,
                                 weights=E(graph)$weight) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (!is.null(weights)) {
    weights <- as.numeric(weights)
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_community_fastgreedy", graph, as.logical(merges),
        as.logical(modularity), weights,
        PACKAGE="igraph0")
} 

community.to.membership <- function(graph, merges, steps, membership=TRUE,
                                    csize=TRUE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  merges <- as.matrix(merges)
  merges <- structure(as.numeric(merges), dim=dim(merges))
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_community_to_membership", graph, merges, as.numeric(steps),
        as.logical(membership), as.logical(csize),
        PACKAGE="igraph0")
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

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  res <- .Call("R_igraph_community_leading_eigenvector_step",
               graph, as.numeric(membership), as.numeric(community),
               options,
               PACKAGE="igraph0")
  class(res) <- "igraph.eigencstep"
  res
}

as.dendrogram.igraph.eigenc <- function(object, hang=-1,
                                         merges=NULL, ...) {
  if (!is.null(merges)) {
    object$merges <- merges
  }
  stopifnot(nrow(object$merges)>0)
  storage.mode(object$merges) <- "integer"
  object$merges <- object$merges + 1L
  if (is.null(object$labels)) {
    labels <- character()
    for (i in 1:max(object$membership+1)) {
      labels[i] <- paste(sep="", "{",
                         paste(collapse=",", which(object$membership==i-1)-1), "}")
    }
    object$labels <- labels
  }
  z <- list()
  object$height <- 1:nrow(object$merges)
  nMerge <- length(oHgt <- object$height)
  if (nMerge != nrow(object$merges))
    stop("'merge' and 'height' do not fit!")
  hMax <- oHgt[nMerge]
  one <- 1:1;
  two <- 2:2 # integer!
  leafs <- nrow(object$merges)+1
  for (k in 1:nMerge) {
    x <- object$merges[k, ]# no sort() anymore!
    if (any(neg <- x < leafs+1))
      h0 <- if (hang < 0) 0 else max(0, oHgt[k] - hang * hMax)
    if (all(neg)) {                  # two leaves
      zk <- as.list(x)
      attr(zk, "members") <- two
      attr(zk, "midpoint") <- 0.5 # mean( c(0,1) )
      objlabels <- object$labels[x]
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
      attr(zk[[2 - isL]], "label") <- object$labels[x[2 - isL]]
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
