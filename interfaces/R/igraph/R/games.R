
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

ba.game <- function(n, power=1, m=NULL, out.dist=NULL, out.seq=NULL,
                    out.pref=FALSE, directed=TRUE, time.window=NULL) {

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
  if (!is.null(time.window) && time.window <= 0) {
    stop("time window size should be positive")
  }
  if (!is.null(m) && m==0) {
    warning("`m' is zero, graph will be empty")
  }
  if (power < 0) {
    warning("`power' is negative")
  }
  
  if (is.null(m) && is.null(out.dist) && is.null(out.seq)) {
    m <- 1
  }
  
  n <- as.numeric(n)
  if (!is.null(m)) { m <- as.numeric(m) }
  if (!is.null(out.dist)) { out.dist <- as.numeric(out.dist) }
  if (!is.null(out.seq)) { out.seq <- as.numeric(out.seq) }
  out.pref <- as.logical(out.pref)

  if (!is.null(out.dist)) {
    out.seq <- as.numeric(sample(0:(length(out.dist)-1), n,
                                 replace=TRUE, prob=out.dist))
  }

  if (is.null(out.seq)) {
    out.seq <- numeric()
  }

  if (!is.null(time.window)) {
    .Call("R_igraph_recent_degree_game", n, power, time.window, m, out.seq,
          out.pref, directed,
          PACKAGE="igraph")
  } else if (power==1) {    
    .Call("R_igraph_barabasi_game", n, m, out.seq, out.pref, directed,
          PACKAGE="igraph")
  } else {
    .Call("R_igraph_nonlinear_barabasi_game", n, power, m, out.seq,
          out.pref, directed,
          PACKAGE="igraph")
  }
}

barabasi.game <- ba.game

erdos.renyi.game <- function(n, p.or.m, type="gnp",
                             directed=FALSE, loops=FALSE, ...) {
  if (is.character(type)) {
    type <- switch(type, "gnp"=0, "gnm"=1)
  }

  .Call("R_igraph_erdos_renyi_game", as.numeric(n), as.numeric(type),
        as.numeric(p.or.m), as.logical(directed), as.logical(loops),
        PACKAGE="igraph")
}

random.graph.game <- erdos.renyi.game

degree.sequence.game <- function(out.deg, in.deg=numeric(0),
                                 method="simple", ...) {

  if (is.character(method)) {
    method <- switch(method, "simple"=0)
  }

  .Call("R_igraph_degree_sequence_game", as.numeric(out.deg),
        as.numeric(in.deg), as.numeric(method),
        PACKAGE="igraph")
}

growing.random.game <- function(n, m=1, directed=TRUE, citation=FALSE) {
  .Call("R_igraph_growing_random_game", as.numeric(n), as.numeric(m),
        as.logical(directed), as.logical(citation),
        PACKAGE="igraph")
}

aging.prefatt.game <- function(graph, n, out.deg=rep(1,n),
                               age.bin=1, pa.exp=1, ac=1,
                               age.exp=1) {

  vc <- vcount(graph)
  ec <- ecount(graph)

  if (vc==0) {
    stop("Cannot start with empty graph")
  }

  born <- as.integer((1:(vc+n)-1) / age.bin)+1
  time <- tail(born, 1)+1
  ind <- c(degree(graph, mode="in"), rep(0, n))
  prob <- (time-born)^-age.exp * (ind+ac)^pa.exp
  edges <- numeric(sum(out.deg)*2)
  edgep <- 1

  for (nn in (vc+1):(vc+n)) {

    sam <- sample(nn-1, out.deg[nn-vc], prob=prob[1:(nn-1)])

    ## aging
    if ( (nn-vc) %% age.bin==0) {
      time <- time + 1
      prob[1:(nn-1)] <- prob[1:(nn-1)] *
        ((time-born[1:(nn-1)]+1)/(time-born[1:(nn-1)]))^-age.exp
    }

    ## increase degree
    for (s in sam) {
      prob[s] <- prob[s] * ((ind[s]+ac+1)/(ind[s]+ac))^pa.exp
      edges[edgep] <- nn ; edgep <- edgep + 1
      edges[edgep] <- s ; edgep <- edgep + 1
    }
  }
  
  res <- add.vertices(graph, n)
  res <- add.edges(res, edges-1)
  res
}

callaway.traits.game <- function(nodes, types, edge.per.step=1,
                                type.dist=rep(1, types),
                                pref.matrix=matrix(1, types, types),
                                directed=FALSE) {

  .Call("R_igraph_callaway_traits_game", as.double(nodes),
        as.double(types), as.double(edge.per.step),
        as.double(type.dist), matrix(as.double(pref.matrix), types, types),
        as.logical(directed))
}

establishment.game <- function(nodes, types, k=1, type.dist=rep(1, types),
                               pref.matrix=matrix(1, types, types),
                               directed=FALSE) {
  .Call("R_igraph_establishment_game", as.double(nodes),
        as.double(types), as.double(k), as.double(type.dist),
        matrix(as.double(pref.matrix), types, types),
        as.logical(directed))
}
