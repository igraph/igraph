
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

membership <- function(communities) {
  if (!is.null(communities$membership)) {
    return(communities$membership)
  } else if (!is.null(communities$merges) &&
             !is.null(communities$modularity)) {
    return(community.to.membership2(communities$merges, communities$vcount,
                                    which.max(communities$modularity)))
  } else {
    stop("Cannot calculate community membership")
  }
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
  }
  invisible(x)
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
  m[el[,1]] != m[el[,2]]
}

is.hierarchical <- function(communities) {
  if (algorithm(communities) %in% c("walktrap", "edge betweenness", "fast greedy", "leading eigenvector")) {
    TRUE
  } else if (algorithm(communities) %in% c("spinglass",
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
    object$labels <- 1:(nrow(object$merges)+1)
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
  one <- 1L
  two <- 2L
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
    res$algorithm  <- "spinglass"
    res$vcount     <- vcount(graph)
    res$membership <- res$membership + 1
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
  res$membership <- res$membership + 1
  res$merges <- res$merges + 1
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
  res$membership <- res$membership + 1
  res$merged <- res$merges + 1
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
  .Call("R_igraph_community_to_membership", graph, merges, as.numeric(steps),
        as.logical(membership), as.logical(csize),
        PACKAGE="igraph")
}

igraph.i.levc.arp <- function(externalP, externalE) {
  f <- function(v) {
    v <- as.numeric(v)
    .Call("R_igraph_i_levc_arp", externalP, externalE, v, PACKAGE="igraph");
  }
  f
}

leading.eigenvector.community <- function(graph, steps=-1, start=NULL,
                                          options=igraph.arpack.default,
                                          callback=NULL, extra=NULL,
                                          env=parent.frame()){

  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  steps <- as.numeric(steps)
  if (!is.null(start)) { start <- as.numeric(start)-1 }
  options.tmp <- igraph.arpack.default; options.tmp[ names(options) ] <- options ; options <- options.tmp
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_community_leading_eigenvector", graph, steps,
               options, start, callback, extra, env,
               environment(igraph.i.levc.arp),
               PACKAGE="igraph")
  res$algorithm <- "leading eigenvector"
  res$vcount <- vcount(graph)
  res$membership <- res$membership + 1
  res$merges <- res$merges + 1
  res$history <- res$history + 1
  class(res) <- "communities"
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
  res$vcount <- vcount(graph)
  res$algorithm <- "multi level"
  res$membership <- res$membership + 1
  res$memberships <- res$memberships + 1
  class(res) <- "communities"
  res
}

optimal.community <- function(graph) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_community_optimal_modularity", graph,
               getIgraphOpt("verbose"), PACKAGE="igraph")
  res$vcount <- vcount(graph)
  res$algorithm <- "optimal"
  res$membership <- res$membership + 1
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

hrg.consensus <- function(graph, hrg=NULL, start=FALSE, num.samples=10000) {
  
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }

  if (is.null(hrg)) {
    hrg <- list(left=c(), right=c(), prob=c(), edges=c(), vertices=c())
  }
  hrg <- lapply(hrg[c("left","right","prob","edges","vertices")], as.numeric)
  start <- as.logical(start)
  num.samples <- as.integer(num.samples)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_hrg_consensus", graph, hrg, start, num.samples,
        PACKAGE="igraph")
  res$parents <- res$parents + 1
  class(res) <- "igraphHRGConsensus"
  class(res$hrg) <- "igraphHRG"
  res
}

hrg.predict <- function(graph, hrg=NULL, start=FALSE, num.samples=10000,
                        num.bins=25) {
  
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  if (is.null(hrg)) { 
    hrg <- list(left=c(), right=c(), prob=c(), edges=c(), 
                vertices=c()) 
  } 
  hrg <- lapply(hrg[c("left","right","prob","edges","vertices")], 
                as.numeric)
  start <- as.logical(start)
  num.samples <- as.integer(num.samples)
  num.bins <- as.integer(num.bins)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_hrg_predict", graph, hrg, start, num.samples,
               num.bins, PACKAGE="igraph")
  res$edges <- matrix(res$edges, ncol=2, byrow=TRUE)
  class(res$hrg) <- "igraphHRG"
  res
}

print.igraphHRG <- function(x, type=c("auto", "tree", "plain"),
                            level=3, ...) {

  type <- igraph.match.arg(type)
  if (type=="auto") {
    type <- if (length(x$left <= 100)) "tree" else "plain"
  }
  if (type=="tree") {
    return(print1.igraphHRG(x, level=level, ...))
  } else {
    return(print2.igraphHRG(x, ...))
  }
}

print1.igraphHRG <- function(x, level=3, ...) {
  cat(sep="", "Hierarchical random graph, at level ", level, ":\n")

  ## Depth of printed top of the dendrogram
  .depth <- function(b, l) {
    l[2] <- max(l[2], nchar(format(x$prob[b], digits=2)))
    if (l[1]==level) { return(l) }
    if (x$left[b] < 0 && x$right[b] < 0) {
      l1 <- .depth(-x$left[b], c(l[1]+1, l[2]))
      l2 <- .depth(-x$right[b], c(l[1]+1, l[2]))
      return(pmax(l1,l2))
    }
    if (x$left[b] < 0)  { return(.depth(-x$left[b],  c(l[1]+1, l[2]))) }
    if (x$right[b] < 0) { return(.depth(-x$right[b], c(l[1]+1, l[2]))) }
    return(l)
  }
  cs <- .depth(1, c(1, 0))
  pw <- cs[2]
  cs <- cs[1] * 4
  vw <- nchar(as.character(length(x$left)+1))
  sp <- paste(collapse="", rep(" ", cs+pw+2+2))

  ## Function to collect all individual vertex children
  .children <- function(b) {
    res <- c()
    if (x$left[b]  < 0) {
      res <- c(res, .children(-x$left[b]))
    } else {
      res <- c(x$left[b]+1, res)
    }
    if (x$right[b] < 0) {
      res <- c(res, .children(-x$right[b]))
    } else {
      res <- c(x$right[b]+1, res)
    }
    return(res)
  }

  ## Recursive printing
  .plot <- function(b, l, ind = "") {
    if (b != 1) {
      he <- format(paste(sep="", ind, "'- B-", b), width=cs)
      ind <- paste("   ", ind)
    } else {
      he <- format(paste(sep="", ind, "B-", b), width=cs)
    }
    ## whether to go left and/or right
    gol <- x$left[b]  < 0 && l < level
    gor <- x$right[b] < 0 && l < level

    ## the children to print
    ch1 <- character()
    if (!gol && x$left[b] < 0) {
      ch1 <- c(ch1, paste(sep="", "B-", -x$left[b]))
    }
    if (!gor && x$right[b] < 0) {
      ch1 <- c(ch1, paste(sep="", "B-", -x$right[b]))
    }
    ch2 <- character()
    if (!gol) {
      if (x$left[b]  <  0) { ch2 <- c(ch2, .children(-x$left[b])) }
      if (x$left[b]  >= 0) { ch2 <- c(ch2, x$left[b] + 1) }
    }
    if (!gor) {
      if (x$right[b] <  0) { ch2 <- c(ch2, .children(-x$right[b])) }
      if (x$right[b] >= 0) { ch2 <- c(ch2, x$right[b] + 1) }
    }

    ## print this line
    lf <- gsub(" ", "x", format(ch2, width=vw), fixed=TRUE)
    lf <- paste(collapse=" ", lf)
    lf <- strwrap(lf, width=getOption("width") - cs - pw - 3 - 2)
    lf <- gsub("x", " ", lf, fixed=TRUE)
    if (length(lf) > 1) {
      lf <- c(lf[1], paste(sp, lf[-1]))
      lf <- paste(collapse="\n", lf)
    }
    op <- paste(sep="", format(he, width=cs),
                " p=", format(x$prob[b], digits=2, width=pw, justify="left"),
                "  ", paste(collapse=" ", lf))
    cat(op, fill=TRUE)

    ## recursive call
    if (x$left[b]  < 0 && l < level) .plot(-x$left[b],  l+1, ind)
    if (x$right[b] < 0 && l < level) .plot(-x$right[b], l+1, ind)
  }

  ## Do it
  if (length(x$left) > 0) .plot(b=1, l=1)

  invisible(x)
}

print2.igraphHRG <- function(x, ...) {
  cat(sep="", "Hierarchical random graph:\n")
  bw <- ceiling(log10(length(x$left)+1))+2
  p <- format(x$prob, digits=1)
  pw <- 4 + max(nchar(p))
  op <- sapply(seq_along(x$left), function(i) {
    lc <- if (x$left[i] < 0) {
      paste(sep="", "B-", -x$left[i])
    } else {
      x$left[i]+1
    }
    rc <- if (x$right[i] < 0) {
      paste(sep="", "B-", -x$right[i])
    } else {
      x$right[i]+1
    }
    paste(sep="", format(paste(sep="", "B-", i), width=bw),
          format(paste(sep="", " p=", p[i]), width=pw),
          "-> ", lc, " ", rc)
  })
  op <- format(op, justify="left")
  cat(op, sep="  ", fill=TRUE)
}

"
## How to print HRGs?

B-1           p=0
'- B-3        p=1  6
   '- B-7     p=1  2
      '- B-5  p=1  1 5
'- B-6        p=1  7
   '- B-2     p=1  4
      '- B-4  p=1  3 8

## The same at levels 1, 2 and 3:

B-1  p=0  B-3 B-6 6 2 1 5 7 4 3 8

B-1     p=0
'+ B-3  p=1  B-7  6 2 1 5
'+ B-6  p=1  B-2  7 4 3 8

B-1        p=0
'- B-3     p=1  6
   '+ B-7  p=1  B-5 2 1 5
'- B-6     p=1  7
   '+ B-2  p=1  B-4 4 3 8

## This can be tedious if the graph is big, as we always have n-1
## internal nodes, we can restrict ourselves to (say) level 3 by default.

## Another possibility is to order the lines according to the group ids.

B-1  p=0  B-3 B-6
B-2  p=1  B-4 4
B-3  p=1  B-7 6
B-4  p=1  3 8
B-5  p=1  1 5
B-6  p=1  B-2 7
B-7  p=1  B-5 2

"
