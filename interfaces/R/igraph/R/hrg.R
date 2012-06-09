
#   IGraph R package
#   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

hrg.fit <- function(graph, hrg=NULL, start=FALSE, steps=0) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  if (is.null(hrg)) { 
    hrg <- list(left=c(), right=c(), prob=c(), edges=c(), 
                vertices=c()) 
  } 
  hrg <- lapply(hrg[c("left","right","prob","edges","vertices")], 
                as.numeric)
  start <- as.logical(start)
  steps <- as.integer(steps)
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_hrg_fit", graph, hrg, start, steps,
               PACKAGE="igraph")
  
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    res$names <- V(graph)$name
  }

  class(res) <- "igraphHRG"
  res
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
  res <- list(consensus=list(parents=res$parents, weights=res$weights),
              hrg=res$hrg)
  class(res$consensus) <- "igraphHRGConsensus"
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

buildMerges <- function(object) {

  ## Build a merge matrix. This is done by a post-order
  ## traversal of the tree.

  S <- numeric()
  vcount <- length(object$left)+1
  nMerge <- vcount-1
  merges <- matrix(0, nrow=vcount-1, ncol=3)
  mptr <- 1
  S[length(S)+1] <- -1
  prev <- NULL
  while (length(S) != 0) {
    curr <- S[length(S)]
    ## coming from parent? going left if possible.
    if (is.null(prev) ||
        (prev < 0 && object$left[-prev] == curr) ||
        (prev < 0 && object$right[-prev] == curr)) {
      if (curr < 0) { S <- c(S, object$left[-curr]) }
    ## coming from left child? going right
    } else if (curr < 0 && object$left[-curr] == prev) {
      S <- c(S, object$right[-curr])
    ## coming from right child? going up
    } else {
      if (curr < 0) {
        merges[mptr,] <- c(object$left[-curr], object$right[-curr], curr)
        mptr <- mptr + 1
      }
      S <- S[-length(S)]
    }
    prev <- curr
  }
  merges
} 

as.dendrogram.igraphHRG <- function(object, hang=0.01, ...) {

  nMerge <- length(object$left)
  merges <- buildMerges(object)
  
  .memberDend <- function(x) {
    r <- attr(x,"x.member")
    if(is.null(r)) {
      r <- attr(x,"members")
      if(is.null(r)) r <- 1:1
    }
    r
  }
  
  oHgt <- 1:nrow(merges)
  hMax <- oHgt[length(oHgt)]
  mynames <- 1:(nMerge+1)
  z <- list()
  
  for (k in 1:nMerge) {
    x <- merges[k,1:2]
    if (any(neg <- x >= 0)) {
      h0 <- if (hang < 0) 0 else max(0, oHgt[k] - hang * hMax)
    }
    if (all(neg)) {                     # two leaves
      zk <- as.list(x+1)
      attr(zk, "members") <- 2L
      attr(zk, "midpoint") <- 1/2       # mean( c(0,1) )
      objlabels <- mynames[x+1]
      attr(zk[[1]], "label") <- objlabels[1]
      attr(zk[[2]], "label") <- objlabels[2]
      attr(zk[[1]], "members") <- attr(zk[[2]], "members") <- 1L
      attr(zk[[1]], "height") <- attr(zk[[2]], "height") <- h0
      attr(zk[[1]], "leaf") <- attr(zk[[2]], "leaf") <- TRUE
    } else if (any(neg)) {              # one leaf, one node
      X <- paste0("g", -x)
      isL <- x[1] >= 0
      zk <- if (isL) list(x[1]+1, z[[X[2]]]) else list(z[[X[1]]], x[2]+1)
      attr(zk, "members") <- attr(z[[X[1+isL]]], "members") + 1L
      attr(zk, "midpoint") <-
        (.memberDend(zk[[1]]) + attr(z[[X[1+isL]]], "midpoint"))/2
      attr(zk[[2 - isL]], "members") <- 1L
      attr(zk[[2 - isL]], "height") <- h0
      attr(zk[[2 - isL]], "label") <- mynames[x[2 - isL]+1]
      attr(zk[[2 - isL]], "leaf") <- TRUE
    } else {                            #two nodes
      X <- paste0("g", -x)
      zk <- list(z[[X[1]]], z[[X[2]]])
      attr(zk, "members") <- attr(z[[X[1]]], "members") +
        attr(z[[X[2]]], "members")
      attr(zk, "midpoint") <- (attr(z[[X[1]]], "members") +
                               attr(z[[X[1]]], "midpoint") +
                               attr(z[[X[2]]], "midpoint"))/2
    }
    attr(zk, "height") <- oHgt[k]
    z[[k <- paste0("g", -merges[k,3])]] <- zk
  }
  z <- z[[k]]
  class(z) <- "dendrogram"
  z
}

as.hclust.igraphHRG <- function(x, ...) {
  merge3 <- buildMerges(x)

  ## We need to rewrite the merge matrix, because hclust assumes
  ## that group ids are assigned in the order of the merges
  map <- order(-merge3[,3])
  
  merge <- merge3[,1:2]
  gs <- which(merge < 0)
  merge[ gs] <- map[ -merge[gs] ]
  merge[-gs] <- -merge[-gs]-1

  ## To get the ordering, we need to recode the merge matrix again,
  ## without using group ids. Here the right node is merged _into_
  ## the left node.
  map2 <- numeric(nrow(merge))
  mergeInto <- merge
  for (i in 1:nrow(merge)) {
    mr <- mergeInto[i,]
    mr[mr > 0] <- -map2[mr[mr>0]]
    mergeInto[i,] <- -mr
    map2[i] <- -mr[1]
  }
  n <- nrow(merge)+1
  hcass <- .Fortran(stats:::C_hcass2, n=as.integer(n),
                    ia=as.integer(mergeInto[,1]),
                    ib=as.integer(mergeInto[,2]),
                    order=integer(n), iia=integer(n), iib=integer(n),
                    PACKAGE="stats")
                    

  mynames <- 1:n
  res <- list(merge=merge, height=1:nrow(merge), order=hcass$order,
              labels=mynames, method=NA_character_,
              dist.method=NA_character_)
  class(res) <- "hclust"
  res  
}

asPhylo.igraphHRG <- function(x, ...) {
  require(ape, quietly=TRUE)
  merge3 <- buildMerges(x)

  ## recode group ids
  map <- order(-merge3[,3])  
  merges <- merge3[,1:2]
  gs <- which(merges < 0)
  merges[ gs] <- map[ -merges[gs] ]
  merges[-gs] <- -merges[-gs]-1

  N <- nrow(merges)
  height <- 1:N
  labels <- 1:(N+1)

  edge <- matrix(0L, 2 * N, 2)
  edge.length <- numeric(2 * N)
  node <- integer(N)
  node[N] <- N + 2L
  cur.nod <- N + 3L
  j <- 1L
  for (i in N:1) {
    edge[j:(j + 1), 1] <- node[i]
    for (l in 1:2) {
      k <- j + l - 1L
      y <- merges[i, l]
      if (y > 0) {
        edge[k, 2] <- node[y] <- cur.nod
        cur.nod <- cur.nod + 1L
        edge.length[k] <- height[i] - height[y]
      } else {
        edge[k, 2] <- -y
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

dendPlot.igraphHRG <- function(x, mode=getIgraphOpt("dend.plot.type"), ...) {

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
    hrgPlotHclust(x, ...)
  } else if (mode=="dendrogram") {
    hrgPlotDendrogram(x, ...)
  } else if (mode=="phylo") {
    hrgPlotPhylo(x, ...)
  }
}

hrgPlotHclust <- function(x, rect=0, colbar=rainbow(rect), hang=.01,
                          ann=FALSE, main="", sub="", xlab="", ylab="",
                          ...) {
  hc <- as.hclust(x)
  ret <- plot(hc, hang=hang, ann=ann, main=main, sub=sub, xlab=xlab,
              ylab=ylab, ...)
  if (rect > 0) {
    rect.hclust(hc, k=rect, border=colbar)
  }
  invisible(ret)
}

hrgPlotDendrogram <- function(x, ...) {
  plot(as.dendrogram(x), ...)
}

hrgPlotPhylo <- function(x, edge.color="#AAAAAAFF", edge.lty=c(1,2), ...) {
  phy <- asPhylo(x)
  plot(phy, edge.color=edge.color, edge.lty=edge.lty[2], ...)
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
  cs <- cs[1] * 3
  vw <- nchar(as.character(length(x$left)+1))
  sp <- paste(collapse="", rep(" ", cs+pw+2+2))
  nn <- if (is.null(x$names)) seq_len(length(x$left)+1) else x$names
  
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
      he <- format(paste(sep="", ind, "'- g", b), width=cs)
      ind <- paste("  ", ind)
    } else {
      he <- format(paste(sep="", ind, "g", b), width=cs)
    }
    ## whether to go left and/or right
    gol <- x$left[b]  < 0 && l < level
    gor <- x$right[b] < 0 && l < level

    ## the children to print
    ch1 <- character()
    if (!gol && x$left[b] < 0) {
      ch1 <- c(ch1, paste(sep="", "g", -x$left[b]))
    }
    if (!gor && x$right[b] < 0) {
      ch1 <- c(ch1, paste(sep="", "g", -x$right[b]))
    }
    ch2 <- numeric()
    if (!gol) {
      if (x$left[b]  <  0) { ch2 <- c(ch2, .children(-x$left[b])) }
      if (x$left[b]  >= 0) { ch2 <- c(ch2, x$left[b] + 1) }
    }
    if (!gor) {
      if (x$right[b] <  0) { ch2 <- c(ch2, .children(-x$right[b])) }
      if (x$right[b] >= 0) { ch2 <- c(ch2, x$right[b] + 1) }
    }

    ## print this line
    ch2 <- as.character(nn[ch2])
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
  cat("Hierarchical random graph:\n")
  bw <- ceiling(log10(length(x$left)+1))+1
  p <- format(x$prob, digits=1)
  pw <- 4 + max(nchar(p))
  nn <- if (is.null(x$names)) seq_len(length(x$left)+1) else x$names
  op <- sapply(seq_along(x$left), function(i) {
    lc <- if (x$left[i] < 0) {
      paste(sep="", "g", -x$left[i])
    } else {
      nn[x$left[i]+1]
    }
    rc <- if (x$right[i] < 0) {
      paste(sep="", "g", -x$right[i])
    } else {
      nn[x$right[i]+1]
    }
    paste(sep="", format(paste(sep="", "g", i), width=bw),
          format(paste(sep="", " p=", p[i]), width=pw),
          "-> ", lc, " ", rc)
  })
  op <- format(op, justify="left")
  cat(op, sep="   ", fill=TRUE)
  invisible(x)
}

# TODO: print as a tree

print.igraphHRGConsensus <- function(x, ...) {
  cat("HRG consensus tree:\n")
  n <- length(x$parents) - length(x$weights)
  id <- c(seq_len(n), paste(sep="", "g", seq_along(x$weights)))
  ch <- tapply(id, x$parents, c)[-1]   # first is zero
  bw <- nchar(as.character(length(x$weights)))
  vw <- max(nchar(id))
  op <- sapply(seq_along(x$weights), function(i) {
    mych <- format(ch[[i]], width=vw)
    if (length(ch[[i]])*(vw+1) + bw + 4 > getOption("width")) {
      mych <- gsub(" ", "x", mych, fixed=TRUE)
      mych <- paste(collapse=" ", mych)
      pref <- paste(collapse="", rep(" ", bw+5))
      mych <- strwrap(mych, width=getOption("width") - bw - 4,
                      initial="", prefix=pref)
      mych <- gsub("x", " ", mych, fixed=TRUE)
      mych <- paste(collapse="\n", mych)
    } else {
      mych <- paste(collapse=" ", mych)
    }
    paste(sep="", "g", format(i, width=bw), " -> ", mych)
  })
  if (max(nchar(op)) < (getOption("width")-4)/2) {
    op <- format(op, justify="left")
    cat(op, sep="   ", fill=TRUE)
  } else {
    cat(op, sep="\n")
  }
  
  invisible(x)
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
