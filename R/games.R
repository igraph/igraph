
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
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
###################################################################

ba.game <- function(n, m=NULL, out.dist=NULL, out.seq=NULL,
                    out.pref=FALSE, ...) {

  # Checks
  if (! is.null(out.seq) && (!is.null(m) || !is.null(out.dist))) {
    warning("if `out.seq' is given `m' and `out.dist' should be NULL")
    m <- out.dist <- NULL
  }
  if (is.null(out.seq) && !is.null(out.dist) && !is.null(m)) {
    warning("if `out.dist' is given `m' will be ignored")
    m <- NULL
  }
  if (!is.null(out.seq) && length(out.seq) != n) {
    stop("`out.seq' should be of length `n'")
  }
  if (!is.null(out.seq) && min(out.seq)<0) {
    stop("negative elements in `out.seq'");
  }
  if (!is.null(m) && m<0) {
    stop("`m' is negative")
  }
  if (!is.null(m) && m==0) {
    warning("`m' is zero, graph will be empty")
  }

  if (is.null(m) && is.null(out.dist) && is.null(out.seq)) {
    m <- 1
  }
  
  n <- as.numeric(n)
  if (!is.null(m)) { m <- as.numeric(m) }
  if (!is.null(out.dist)) { out.dist <- as.numeric(out.dist) }
  if (!is.null(out.dist)) { out.seq <- as.numeric(out.seq) }
  out.pref <- as.logical(out.pref)

  if (!is.null(out.dist)) {
    out.seq <- as.numeric(sample(0:(length(out.dist)-1), n,
                                 replace=TRUE, prob=out.dist))
  }
  
  res <- graph.empty(n=n, ...)
  res <- add.edges(res, .Call("REST_ba_game", n, m, out.seq, out.pref,
                              PACKAGE="igraph"))
  
  res
}

barabasi.game <- ba.game

erdos.renyi.game <- function(n, p, directed=FALSE, loops=FALSE, ...) {

  n <- as.numeric(n)
  p <- as.numeric(p)

  if (n<0) {
    stop("`n' should be positive")
  }
  if (p<0 || p>1) {
    stop("`p' should be between zero and one")
  }

  if (p==0 || n==1) {
    res <- graph.empty(n=n, directed=directed, ...)
  } else {
    # generate enough waiting times
    possible.edges <-
      if ( directed && loops) n**2
      else if ( directed && !loops) n*(n-1)
      else if (!directed && loops) n*(n+1)/2
      else n*(n-1)/2
    s <- cumsum(rgeom(possible.edges*p*1.1, p)+1)-1
    while (s[length(s)] < possible.edges) {
      more <- cumsum(rgeom(possible.edges*p*0.05, p)+1)+s[length(s)]
      s <- c(s, more)
    }
    s <- s [ s <= possible.edges ]

    # ok, calculate the edges
    if (directed && loops) {
      from <- s %/% n + 1
      to <- s %% n + 1
    } else if (directed && !loops) {
      from <- s %/% n + 1
      to <- s %% (n-1) + 1
      to [ from == to ] <- n
    } else if (!directed && loops) {
      from <- ceiling((sqrt(8*s+1)-1)/2)
      to <- s-from*(from-1)/2
    } else {
      from <- ceiling((sqrt(8*s+1)-1)/2)+1
      to <- s-(from-1)*(from-2)/2
    }
    
    # ok, create graph
    res <- graph( as.numeric(t(matrix(c(from, to), nc=2))),
                 directed=directed, n=n, ...)
  }
  
  res
}

random.graph.game <- erdos.renyi.game

degree.sequence.game <- function(out.deg, in.deg=NULL, method="simple", ...) {

  if (any(out.deg < 0)) {
    stop("Negative degree not allowed in `out.deg'")
  }
  if (!is.null(in.deg) && any(in.deg<0)) {
    stop("Negative degree not allowed in `in.deg'")
  }
  if (is.null(in.deg) && sum(out.deg) %% 2 != 0) {
    stop("Total degree should be even")
  }
  if (!is.null(in.deg) && length(out.deg) != length(in.deg)) {
    stop("Length of `in.deg' should match length of `out.deg'")
  }
  if (!is.null(in.deg) && sum(in.deg) != sum(out.deg)) {
    stop("Total in-degree should match total out-degree")
  }
  if (method != "simple") {
    stop("Invalid `method', see docs")
  }

  if (is.null(in.deg)) {
    directed <- FALSE
    edges <- rep(1:length(out.deg), times=out.deg)
    edges <- sample(edges, length(edges))
  } else {
    directed <- TRUE
    from <- rep(1:length(out.deg), times=out.deg)
    to <- rep(1:length(in.deg), times=in.deg)
    from <- sample(from, length(from))
    to <- sample(to, length(to))
    edges <- as.numeric(t(matrix(c(from, to), nc=2)))
  }

  res <- graph.empty(n=length(out.deg), directed=directed, ...)
  res <- add.edges(res, edges)  
  
  res
}
