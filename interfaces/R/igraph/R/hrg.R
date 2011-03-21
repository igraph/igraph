
#   IGraph R package
#   Copyright (C) 2011  Gabor Csardi <csardi.gabor@gmail.com>
#   Rue de l'Industrie 5, 1005 Lausanne, Switzerland
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
  cat("Hierarchical random graph:\n")
  bw <- ceiling(log10(length(x$left)+1))+1
  p <- format(x$prob, digits=1)
  pw <- 4 + max(nchar(p))
  op <- sapply(seq_along(x$left), function(i) {
    lc <- if (x$left[i] < 0) {
      paste(sep="", "g", -x$left[i])
    } else {
      x$left[i]+1
    }
    rc <- if (x$right[i] < 0) {
      paste(sep="", "g", -x$right[i])
    } else {
      x$right[i]+1
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
